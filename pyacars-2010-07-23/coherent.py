#!/usr/bin/env python
# -*- coding: utf-8 -*-

from gnuradio import gr
from gnuradio import audio
from gnuradio import trellis
import os
from sys import stdout, stderr, exit
from cmath import rect, pi
from math import log, sin, cos, log10
from acars_message import   message, frames_recoverer, \
                            hr_output_formater, csv_output_formater, \
                            raw_output_formater, \
                            MessageParseError
from ambiguities_resolver import ambiguities_resolver
from threading import Thread
from optparse import OptionParser
from datetime import datetime
from locale import setlocale, LC_ALL

DEBUG=False

if DEBUG:
    import rpdb2; rpdb2.start_embedded_debugger("123")

def sinc(x):
    if x == 0:
        return 1
    return sin(pi * x) / (pi * x)


class acars(gr.hier_block2):
    def __init__(self,input_rate):
        
        gr.hier_block2.__init__(self, "Coherent ACARS demodulator",
                gr.io_signature(1, 1, gr.sizeof_float),
                gr.io_signature(1, 1, gr.sizeof_char))
        
        self.baudratre = 2400.
        self.excess_smpl = 1.7

        # Input filter design params.
        self.input_rate = float(input_rate)
        self.infilt_cutoff_low = 300
        self.infilt_cutoff_high = 3300
        self.infilt_trans_width = 200
        self.infilt_att = 60

        #
        self.cma_numtaps = 1
        self.cma_mu = .06           # The constatnt modulus adaptive equalizer 
                                    # mu-constant.
                                    # MonteCarlo
        
        # Costas loops params.
        self.costas1_alpha = 0.05   # 1st order gain, typically 0.01-0.2
        self.costas1_beta = self.costas1_alpha**2/16  
                                    # 2nd order gain, typically alpha^2/4.0
        self.costas2_alpha = 0.05   # 1st order gain, typically 0.01-0.2
        self.costas2_beta = self.costas1_alpha**2/16  
                                    # 2nd order gain, typically alpha^2/4.0
        self.costas_order = 2
        
        self.pfb_alpha = 12         # Polyphase filterbank clock recovery alpha.

        self.pfb_interp = 64        # Polyphase filterbank clock recovery 
                                    # interpolation factor.
                                    
        self.pfb_symbs = 9          # How much of the symbols the interpolation 
                                    # filter interpolates over -- the length 
                                    # (in symbols) of the impulse response of 
                                    # the pfb filters
                                     
        self.pfb_max_rate_deviation = 0.05
                                    # Too big value seems to cause 
                                    # "hangs" on some messages. PFB.
                                    # In such hangs, pfb outputs 
                                    # periodically, but the samples 
                                    # are not recognized. Maybe 
                                    # they are all same --- zeros 
                                    # probably.

        self.vitter_k = 2**8        # The count of symbols which the Vitter 
                                    # algorithm will be computing the path costs
                                    # from.

        self.mafilt_recursive = 1   # Use recursive MA filters.

    
    def mk_debug_sink(self, block, name, sink):
        self.connect(block, sink)
        setattr(self, name, sink)


    def debugging_connects(self):
        self._PS_sink = gr.vector_sink_c()
        self._AGC_sink = gr.vector_sink_c()
        self._REC_sink = gr.vector_sink_c()
        self._HILB_sink = gr.vector_sink_c()
        self._PLL1_REF_sink = gr.vector_sink_c()
        self._PLL2_REF_sink = gr.vector_sink_c()
        self._PLL1_OUT_sink = gr.vector_sink_c()
        self._PLL2_OUT_sink = gr.vector_sink_c()
        self._PFB_ERR_sink = gr.vector_sink_f()
        self._PFB_RATE_sink= gr.vector_sink_f()
        self._PFB_K_sink= gr.vector_sink_f()


        tb.connect(self._psf, self._PS_sink)
        tb.connect(self._agc, self._AGC_sink)
        tb.connect(self._clock_recovery, self._REC_sink)
        tb.connect(self._hilbfilt, self._HILB_sink)
        tb.connect((self._pll1,1), self._PLL1_REF_sink)
        tb.connect((self._pll2,1), self._PLL2_REF_sink)
        tb.connect((self._pll1,0), self._PLL1_OUT_sink)
        tb.connect((self._pll2,0), self._PLL2_OUT_sink)
        tb.connect((self._clock_recovery,1), self._PFB_ERR_sink)
        tb.connect((self._clock_recovery,2), self._PFB_RATE_sink)
        tb.connect((self._clock_recovery,3), self._PFB_K_sink)


    def build_graph (self):
        # The nyquinst kriterium.

        if self.infilt_trans_width <=0:
            raise ValueError("infilt_trans_width must be positive")
        if self.infilt_cutoff_high < self.infilt_cutoff_low:
            raise ValueError("infilt_cutoff_low must be lesser "
            "than infilt_cutoff_high")
        if self.costas1_alpha<=0:
            raise ValueError("costas1_alpha must be positive")
        if self.costas2_alpha<=0:
            raise ValueError("costas2_alpha must be positive")
        if self.vitter_k <= 0:
            raise ValueError("vitter_k must be positive")

        infilt_maxfreq = self.infilt_cutoff_high + self.infilt_trans_width
        if infilt_maxfreq * 2 >= self.input_rate:
            raise Error("Sample_rate too low.")

        self._input_coefs = input_coefs = \
            gr.firdes.complex_band_pass_2 (
                    1, # gain
                    self.input_rate,
                    self.infilt_cutoff_low,
                    self.infilt_cutoff_high,
                    self.infilt_trans_width,
                    self.infilt_att )
        
        if self.excess_smpl<1:
            raise ValueError("excess_smpl cannot be < 1")
        desired_minfreq = self.infilt_cutoff_high * self.excess_smpl
        self.decimation = decimation = int( self.input_rate // desired_minfreq )

        if decimation > 2:
            print >>stderr, "Input will be decimated by factor %d. " \
                    "You may want to save some power --- lower the input " \
                    "sampling rate a bit." % (decimation)
        else:
            print >>stderr, "Using decimation factor %d." % decimation

        self._hilbfilt = hilbfilt = gr.fir_filter_fcc (decimation, input_coefs)
        
        self._agc = agc = gr.cma_equalizer_cc(
                self.cma_numtaps, 1, self.cma_mu    )

        pllhalfbw = 1200 * 200e-6   # By Arnic
        pllhalfbw += 3
        
        self._pll1 = pll1 = gr.costas_loop_cc(
                self.costas1_alpha,
                self.costas1_beta,
                max_freq=(1200+pllhalfbw)*2*pi/self.input_rate*decimation,
                min_freq=(1200-pllhalfbw)*2*pi/self.input_rate*decimation,
                order = self.costas_order)
        self._pll2 = pll2 = gr.costas_loop_cc(
                self.costas2_alpha,
                self.costas2_beta,
                max_freq=(2400+pllhalfbw*2)*2*pi/self.input_rate*decimation,
                min_freq=(2400-pllhalfbw*2)*2*pi/self.input_rate*decimation,
                order = self.costas_order)

        
        T = self.input_rate/self.baudratre/self.decimation # samples per symbol
        if 1 == self.mafilt_recursive:
            self._psf = gr.moving_average_cc( int(T*2), 1./T/2 )
        elif 0 == self.mafilt_recursive:
            self._psf = gr.fir_filter_ccf(1, [1./T/2] * int(T*2))
        else:
            raise ValueError("mafilt_recursive must be 1 or 0.")
        

        psf = self._psf

        alpha = self.pfb_alpha 
        interp = self.pfb_interp        
        symbs = self.pfb_symbs          # How much symbols convolve over.
        ntaps = interp*symbs*T          # How much taps of prototype filter
        sinc_time = range(-int(ntaps//2), int(ntaps//2))
        # Sinc should go trought zeros on symbol sampling times, 
        # except of the currently sampled symbol.
        norm = 2. / T
        taps = [ sinc(2.*x * symbs / len(sinc_time)) / norm \
                        for x in sinc_time ]
        self._clock_recovery = clock_recovery \
                = gr.pfb_clock_sync_ccf(
                        T,
                        alpha, 
                        taps=tuple(taps), 
                        filter_size=interp,
                        max_rate_deviation=self.pfb_max_rate_deviation
                        )

        if DEBUG:
            a=[]
            for i in range(interp):
                a+= self._clock_recovery.channel_taps(i),
            print len(a)
            from pylab import imshow, show, array, \
                              subplot, plot, fft, unwrap, angle
            subplot(121)
            imshow(array(a), interpolation="nearest")
            subplot(122)
            plot(array(a).transpose())
            show()


        
#        Surely, the psf could be removed. The PFB could do all the work, but it 
#        doesn't work now. I can't see the reason.

#        psf = self._psf = gr.kludge_copy(gr.sizeof_gr_complex)
#        ntaps = 4 * T * interp
#        taps = [0]              * int(T*2*interp) \
#             + [1./T/2/interp]  * int(T*2*interp)
#        print taps, len(taps), ntaps

        
        f2c = gr.float_to_complex()

        self.connect (self,  hilbfilt, agc)
        self.connect (agc, pll1, gr.complex_to_real(), (f2c, 1))
        self.connect (agc, pll2, gr.complex_to_real(), (f2c, 0))
        if self.mafilt_recursive:
            MA_MAXITER = 200
        else:
            MA_MAXITER = 0

        self._ambres = ambres = ambiguities_resolver(DEBUG=False, MA_MAXITER=0)
        self.connect (
                f2c,
                psf,
                clock_recovery,
                ambres,
        )
        
        try:
            fsm = trellis.fsm(os.environ.get("PYACARS_PATH", '') + "fsm")
        except Exception as e: # Just for any case...
            sys.exit(3)

        vit_params = (
                fsm,
                self.vitter_k,
                -1,     # Start state.
                -1,     # Final state.
                1,      # Dimensionality.
                [ -1, -1j, +1j, 1 ],    # Constellation
                trellis.TRELLIS_EUCLIDEAN
#                trellis.TRELLIS_HARD_SYMBOL
#                trellis.TRELLIS_HARD_BIT
#                trellis.TRELLIS_MIN_SUM
#                trellis.TRELLIS_SUM_PRODUCT
        )
        vitc = trellis.viterbi_combined_cb(*vit_params)
        self._vitc = vitc
        
        #                         SYN      SOH    
        # When start_d is 1, the "01101000 1000000" access key works good and 
        # produces 171 CRC OK messages and 2000 parity error messages.
        #
        # When start_d is 0, the "00 1000000" or "0 1000000" works good.
        # produces 177 CRC OK messages and 17000 parity error messages.
        # 
        # Because the number of parity error messages doesn't take a big part
        # in total demodulation time, it may be sacrfised for the number of 
        # correctly demodulated messages. 
        #
        # I'm unsure about which access code is the best.
        #
        startcode = (
                ""         # last bits of SYN
                "01000000"   # SOH without the last bit. The last one will be 
                            # checked by message framer.
                )
        start_d = 0         # If max difference is too big, we may start too 
                            # soon. 1 looks like a good value.

        self.connect(
                ambres,
                (vitc,0),
                gr.correlate_access_code_bb(startcode, start_d),
                self
        )


# TODO Deadlock may occur if itering thread dies.
def file_ds_glue(file_ds, buff_size):
    while True:
        buf = os.read(file_ds, buff_size)
        if buf == "":
            break
        for c in buf:
            yield ord(c)


class messages_displayer(Thread):
    def __init__(self, messages_iter, summary_pos, 
                       summary_file, formatter):
        Thread.__init__(self)
        self.msg_it = messages_iter
        self.summary_pos = summary_pos
        self.summary_file = summary_file
        self.formatter = formatter


    def display_summary(self, stream):
        fr = self._fr
        print >>stream, "Summary:"
        print >>stream, "%6d messages with both parity and CRC correct." \
                % fr.correct_counter
        print >>stream, "%6d messages with bad parity, CRC not checked." \
                % fr.parity_errs_counter
        print >>stream, "%6d messages with good parity but bad CRC." \
                % fr.crc_miss_counter


    def run(self):
        correct_counter = crc_miss_counter = parity_errs_counter = 0
        self._fr = fr = frames_recoverer(self.msg_it)
        self.formatter.write_header()
        for orig in fr:
            try:
                self.orig = orig
                m = message(orig, datetime.utcnow())
            except MessageParseError as e:
                print >>stderr, e
                continue
            self.formatter.write(m)
            if self.summary_pos == 'msg':
                self.display_summary(self.summary_file)
                print
        self.formatter.write_footer()
        if self.summary_pos == 'end':
            self.display_summary(self.summary_file)


if __name__ == '__main__':
    # Let the python get locale.
    # It is needed to get the localized date in human-readable messages.
    setlocale(LC_ALL, '')

    # build the flow graph
    tb = gr.top_block()

    parser = OptionParser()
    
    usage = "usage: %prog { wav WAV_FILE | stdin SAMPLE_RATE | " \
            "rec SAMPLE_RATE [ DEVICE_NAME ] }"
    parser = OptionParser(usage=usage)

    parser.set_defaults( summary_pos="end", summary_file="", out_form="hr" )
    
    parser.add_option(
            "-f", dest="out_form",
            help =  "The format of the output. "
                    "csv - comma-separated values (Excel like dialect), "
                    "hr - human readable format, "
                    "raw - decoding date in ISO8601 format, followed by the raw message.")
    
    parser.add_option(
            "-s", dest="summary_pos",
            help="Where to print the summary. "
            "'msg' for after message, 'end' for the end"   )
    
    parser.add_option(
            "-S", dest="summary_file",
            help="File to open and write the summary to. "
            "Empty string means a stderr."
    )
    
    parser.add_option(
            "--set_status", dest="setstatus",
            action="store_true"
    )
    
    dem_options = {
            'excess_smpl' : 'float',
            'cma_numtaps' : 'int',
            'cma_mu' : 'float',
            'pfb_alpha': 'float',
            'pfb_interp': 'int',
            'pfb_symbs': 'int',
            'vitter_k': 'int',
            'infilt_cutoff_low': 'float',
            'infilt_cutoff_high': 'float',
            'infilt_trans_width': 'float',
            'costas1_alpha': 'float',
            'costas2_alpha': 'float',
            'costas_order': 'int',
            'mafilt_recursive': 'int',
    }
    for optname, type in dem_options.iteritems():
        parser.add_option("--" + optname, dest=optname, type=type )
    
    (options, args) = parser.parse_args()

    if not len(args) > 0:
        parser.error("Please, specify the INPUT_TYPE.")
    input_type = args.pop(0)

    if input_type == "wav":
        if not len(args) > 0:
            parser.error("Please, specify the WAV_FILE.")
        src = gr.wavfile_source(args.pop(0), False)
        samplerate = src.sample_rate()
    elif input_type == "stdin":
        if not len(args) > 0:
            parser.error("Please, specify the SAMPLE_RATE.")
        samplerate=int(args.pop(0))
        src = gr.file_descriptor_source(4, 0, False)
    elif input_type == "rec":
        if not len(args) > 0:
            parser.error("Please, specify the SAMPLE_RATE.")
        samplerate=int(args.pop(0))
        if not len(args) > 0:
            devname='default'
        else:
            devname = args.pop(0)
        src = audio.source(samplerate, devname)
    else:
        parser.error("INPUT_TYPE must be one of: wav, stdin, rec.")

    if len(args) != 0:
        print >> stderr, "Warning: some extra arguments detected. Will be ignored."

    if samplerate < 3000*2:
        print >> stderr, "Sample rate is too low."

    r_fd, w_fd = os.pipe()
    f = acars(samplerate)

    for optname in dem_options:
        val = getattr(options, optname)
        if val is not None:
            setattr(f, optname, val)
    f.build_graph()
    tb.connect(src, f, gr.file_descriptor_sink(gr.sizeof_char, w_fd))
    
    DEBUG = False
    if DEBUG:
        f.debugging_connects()

    tb.start()


    # Prepare the interface between the messages_displayer and the 
    # reading end of the pipe from the GNU Radio. 
    _file_ds_glue = file_ds_glue(r_fd, 2**10)
    
    # Prepare the selected output formatter.
    FORMATTERS_MAP = {  "csv" : csv_output_formater,
                        "hr"  : hr_output_formater,
                        "raw" : raw_output_formater }
    try:
        formatter = FORMATTERS_MAP[options.out_form](stdout)
    except KeyError:
        raise Exception("Unknown output formatter specified!")
    
    # Prepare the output stream for the summary to be printed to.
    summary_file_name = options.summary_file
    if "" == summary_file_name:
        summary_file = stderr
    else:
        summary_file = open(summary_file_name, "w")

    # Make the messages_displayer with prepared parameters and start it.
    msgdsp = messages_displayer(_file_ds_glue, options.summary_pos, 
            summary_file, formatter)
    msgdsp.start()
    
    # Wait for processing until EOF.
    tb.wait()
    os.close(w_fd)  #  w_fd must be closed when done to signal the r_fd reader that 
                    # it is done as well.

    msgdsp.join()   # Wait for the demodulated but not displayed messages in 
                    # the pipe.
    os.close(r_fd)
    
    if DEBUG:
        from numpy import array
        PS = array(f._PS_sink.data())
        HILB = array(f._HILB_sink.data())
        PLL1_OUT = array(f._PLL1_OUT_sink.data())
        PLL2_OUT = array(f._PLL2_OUT_sink.data())
        PLL1_REF = array(f._PLL1_REF_sink.data())
        PLL2_REF = array(f._PLL2_REF_sink.data())
        AGC = array(f._AGC_sink.data())
        REC = array(f._REC_sink.data())
        PFB_ERR = array(f._PFB_ERR_sink.data())
        PFB_RATE = array(f._PFB_RATE_sink.data())
        PFB_K = array(f._PFB_K_sink.data())

    # Close the summary stream. Note the summary_file may be the stderr!
#    summary_file.close()
    
    if options.setstatus:
        exit(msgdsp._fr.correct_counter)
