#!/usr/bin/env python
# -*- coding: utf-8 -*-

from gnuradio import gr
from cmath import rect, pi


class ambiguities_resolver( gr.hier_block2 ):
    """
    Resolves Costas loops ambiguities using knowledge of ACARS message preamble.

    The leading part of ACARS message consits of sequence of ones folowed by 
    0101010101010 bit seqence, which will be transmitted as high frequency 
    followed by a low frequency. In constellation diagram, the high frequency is 
    ±f_h followed by ±f_l. When using two costas loops to detect the frequencies 
    present, their phase is ambiguous, but using correlation (moving averge), 
    it can be recovered.

    Input is stream of symbols.

    Output 0 is the correctly flipped symbols. Flip takes place just before the 
    pre-key ends and the 0101010101010 sequence starts.

    Output 1 is the absolute value of pre-key correlator output. Usable to 
    detect a message on the channel.

    Output 2 is one when the 0101010101010 sequence was found just after the 
    pre-key thus the output will be flipped correctly.
    """

    def __init__(self, len_fh=70, len_fl=13, MA_MAXITER=4096, DEBUG=False):

        gr.hier_block2.__init__(self, "ACARS MSK amabiguities resolver",
                gr.io_signature(1, 1, gr.sizeof_gr_complex),
                gr.io_signature(3, 3, gr.sizeof_gr_complex))

        
        # Prepare blocks.

        in2R, in2I = gr.complex_to_real(), gr.complex_to_imag()
        # No symbols are actually needed. We just need to slice them.
        symout = [ 0, 1, 2, 3 ] 
        constell = [ rect(1, pi/2 * i) for i in symout ]
        constell_dec = gr.constellation_decoder_cb(constell, symout)
        symbs = gr.chunks_to_symbols_bc(constell)
        # We will need to separate I and Q branches.
        sym2r, sym2i = gr.complex_to_real(), gr.complex_to_imag()
        # Pattern correlators -- moving averges.
        if 0 < MA_MAXITER:
            ma_R = gr.moving_average_ff(len_fh, 2./(len_fh+len_fl), MA_MAXITER)
            ma_I = gr.moving_average_ff(len_fl, 2./(len_fh+len_fl), MA_MAXITER)
        else:
            ma_R = gr.fir_filter_fff(1, [2./(len_fh+len_fl)]*len_fh)
            ma_I = gr.fir_filter_fff(1, [2./(len_fh+len_fl)]*len_fl)
        # We will need to compensate for ma_I delay.
        in2R_delay = gr.delay(gr.sizeof_float, len_fl)
        in2I_delay = gr.delay(gr.sizeof_float, len_fl)
        corr_delay = gr.delay(gr.sizeof_float, len_fl)
        th_out_R, th_out_I = gr.threshold_ff(-0,0), gr.threshold_ff(-0,0)
        snh_R, snh_I = gr.sample_and_hold_ff(), gr.sample_and_hold_ff()
        # Symbol flippers.
        mul_out_R, mul_out_I = gr.multiply_ff(), gr.multiply_ff()
        C = gr.float_to_complex()

        # For computing the absolute value.
        mul_corr_R, mul_corr_I = gr.multiply_ff(), gr.multiply_ff()
        
        # Prekey correlating creates a hill, a 101010... correlator creates 
        # a peak on that hill when summed together to make the peak visible 
        # over the noise.
        sum = gr.add_ff()
        peak = gr.peak_detector_fb(
                threshold_factor_rise=1,
                threshold_factor_fall=1,
                look_ahead=30,
                alpha=.1
        )
        floor_th = gr.threshold_ff(1,1)
        slicer = gr.binary_slicer_fb()
        and_block = gr.and_bb()

        # Wiring.
        self.connect(
                self,
                constell_dec,
                symbs,
                sym2r,
                ma_R,
                corr_delay,
                th_out_R,
                # {0,1} -> {-1,1}
                gr.multiply_const_ff(2),
                gr.add_const_ff(-1),
                snh_R,
                mul_out_R,
                (C,0)
        )
        self.connect(
                self,
                in2R,
                in2R_delay,
                (mul_out_R,1)
        )
        self.connect(
                symbs,
                sym2i,
                ma_I,
                th_out_I,
                # {0,1} -> {-1,1}
                gr.multiply_const_ff(2),
                gr.add_const_ff(-1),
                snh_I,
                mul_out_I,
                (C,1)
        )
        self.connect(
                self,
                in2I,
                in2I_delay,
                (mul_out_I,1)
        )
        self.connect(
                corr_delay,
                mul_corr_R,
                sum,
                peak,
                and_block,
                (snh_R,1)
        )
        # Absolute value computing hack.
        self.connect(
                corr_delay,
                gr.threshold_ff(0,0), 
                # {0,1} -> {-1,1}
                gr.multiply_const_ff(2),
                gr.add_const_ff(-1),

                (mul_corr_R,1)
        )
        self.connect(
                ma_I,
                mul_corr_I,
                (sum,1)
        )
        # Absolute value computing hack.
        self.connect(
                ma_I, 
                gr.threshold_ff(0,0), 
                # {0,1} -> {-1,1}
                gr.multiply_const_ff(2),
                gr.add_const_ff(-1),
                (mul_corr_I,1)
        )
        self.connect(
                sum,
                floor_th,
                # {0,1} -> {-.5,.5}, because slicer(0) == 1
                gr.add_const_ff(-.5),
                slicer,
                (and_block,1)
        )
        self.connect(
                and_block,
                (snh_I,1)
        )
        
        self.connect(C, gr.conjugate_cc(), self)
        self.connect(ma_I, (self,1))
        self.connect(and_block, (self,2))
        
        if DEBUG:
            self.mk_debug_sink (symbs, "SYMBS", gr.vector_sink_c())
            self.mk_debug_sink (sym2r, "S2R", gr.vector_sink_f())
            self.mk_debug_sink (sym2i, "S2I", gr.vector_sink_f())
            self.mk_debug_sink (peak, "PEAK", gr.vector_sink_b())
            self.mk_debug_sink (and_block, "AND", gr.vector_sink_b())
            self.mk_debug_sink (C, "C", gr.vector_sink_c())
            self.mk_debug_sink (corr_delay, "D", gr.vector_sink_f())
            self.mk_debug_sink (ma_R, "MA_R", gr.vector_sink_f())
            self.mk_debug_sink (ma_I, "MA_I", gr.vector_sink_f())
            self.mk_debug_sink (snh_I, "SNH_I", gr.vector_sink_f())
            self.mk_debug_sink (snh_R, "SNH_R", gr.vector_sink_f())
            self.mk_debug_sink (sum, "SUM", gr.vector_sink_f())

    
    def mk_debug_sink(self, block, name, sink):
        self.connect(block, sink)
        setattr(self, name, sink)
