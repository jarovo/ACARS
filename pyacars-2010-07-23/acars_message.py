#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#   Copyleft @ Jaroslav Henner.
#   All wrongs reserved.
#

import re
from csv import reader, QUOTE_ALL, DictWriter
from copy import copy
from sys import stderr, exit, exit
import os



class MessageParseError( Exception ):
    def __init__( self, raw_message ):
        self.raw_message = raw_message
        Exception.__init__( self, 
                "Message parse error. Raw messge: %r" % raw_message )
    
    
def parity( int_type ):
    """
    From http://wiki.python.org/moin/BitManipulation. 

    Returns 0 if there are an even number of set bits, and -1 if there are an 
    odd number.

    Time complexity = O(n), where n is a number of ones in the int_type.
    """
    parity = 0
    while int_type:
        parity = ~parity
        int_type = int_type & (int_type - 1)
    return( parity )
    

def receive_byte_lsb_first( message_iter, init_val, bits ):
    """
    Collects LSB from values readen using the message_iter. Bits are shifted to 
    the output byte from left to right. The output byte is initialized to 
    init_val prior any shifting.

    @param message_iter Iterator to the message stream for reading the bits from.
    @param init_val     Initial value for the output byte
    @param bits         Number of bits to be read from the message_iter.
    """
    byte = init_val & 0xff
    for i in range( bits ):
        byte >>= 1
        byte |= ( next( message_iter ) << 7 ) & 0x80
    return byte


def to_char( byte ):
    byte &= 0xff >> 1
    return byte


MESSAGE_OK = 'Message OK'
CRC_MISMATCH = 'CRC mismatch!'
PARITY_MISMATCH = 'Parity mismatch!'


def update_CRC_byte( byte, crc ):
    crc ^= byte
    for i in range(8):
        if crc & 0x0001:
            crc >>= 1
            crc ^= 0x8408 # crc ^= 0x1021
        else:
            crc >>= 1
    return crc


def decode_message( data_it ):
    message = [ 1 ]     # Fake the SOH. It actually is recieved, but not by 
                        # decode_message().
    computed_crc = 0

    while True:
        b = receive_byte_lsb_first( data_it, 0, 8 )
        c = to_char( b )
        message += ( c, )
        computed_crc = update_CRC_byte(b, computed_crc)            
        ETX = 0b00000011
        ETB = 0b00010111
        if  ETX == c or ETB == c:
            message = "".join(map(chr,message))

            received_crc = 0
            received_crc |= receive_byte_lsb_first( data_it, 0, 8 ) << 0
            received_crc |= receive_byte_lsb_first( data_it, 0, 8 ) << 8

#            if 0x7f != receive_byte_lsb_first( data_it, 0, 8 ):
#                print "Not a properly terminated message!"
#                print message
                # TODO Do we have to take care?

#            print "received %x, %x computed" % (received_crc, computed_crc)
            if received_crc != computed_crc:
#                print " CRC miss!"
                return message, CRC_MISMATCH, data_it
            else:
#                print " CRC OK."
                return message, MESSAGE_OK, data_it

        if -1 != parity( b ):
            message = "".join(map(chr,message))
            return message, PARITY_MISMATCH, data_it


class frames_recoverer:
    def __init__( self, channel_data, bad_parity_callback=None ):
        self.channel_data = channel_data
        self.channel_data_it = iter( channel_data )
        self.correct_counter = 0
        self.parity_errs_counter = 0
        self.crc_miss_counter = 0
        self.crc = 0
        self.bad_parity_callback = bad_parity_callback

    
    def __iter__( self ):
        return self


    def find_start( self ):
        channel_data_it = self.channel_data_it
        while True:
            # The last bit of SOH in network order is 0 and it is marked ...
            if 0b10 == next( channel_data_it ):
                # ... iterator now is on the byte right after SOH.
                return self.channel_data_it


    def next( self ):
        while True:
            self.find_start()
            message, correctness, end_it \
                    = decode_message( self.channel_data_it )
            if correctness is MESSAGE_OK:
                self.channel_data_it = end_it
                self.correct_counter += 1
                return message
            elif correctness is PARITY_MISMATCH:
                self.parity_errs_counter += 1
                if self.bad_parity_callback is not None:
                    self.bad_parity_callback(message)
            elif correctness is CRC_MISMATCH:
                self.crc_miss_counter += 1
            else:
                assert(False)


def prepare_label_dict(label_dict_fname):
    """
    Returns a labels dictionary readen from csv file.
    """
    with open(label_dict_fname, "rb") as ld:
        r = reader( ld, 
                    delimiter=":", 
                    quotechar="'", 
                    quoting=QUOTE_ALL, 
                    skipinitialspace=True)
        l = ( (key, (direction, text)) for key, direction, text in r )
        return dict(l)


def prepare_originator_dict(originator_dict_fname):
    """
    Returns originator dictionary readen from csv file.
    """
    with open(originator_dict_fname, "rb") as od:
        r = reader( od, 
                    delimiter=":", 
                    quotechar="'", 
                    quoting=QUOTE_ALL, 
                    skipinitialspace=True)
        return dict(r)


# TODO Pack into the object everything to remove this initialisation, which 
# cannot be parametrized.
try:
    _LABEL_DICT = prepare_label_dict(
            os.environ.get("PYACARS_PATH", '') + "label.dict"   )
    _ORIGINATOR_DICT = prepare_originator_dict(
            os.environ.get("PYACARS_PATH", '') + "originator.dict"    )
except IOError as e:
    print >>stderr, e
    print >>stderr, "Is the PYACARS_PATH environment variable set " \
                    "properly (eg.:/usr/local/pyacars/)?"
    exit(2)
    
def decomposed_label( label ):
    return _LABEL_DICT.get(label, ('', label))



_REGEXP = '^[\x01](?P<MODE>.)(?P<ADDRESS>.{7})(?P<TECH_ACK>.)' \
        + '(?P<LABEL>..)(?P<UBI_DBI>.)?(?:\x02(?P<TEXT>.*))?' \
        + '(?P<END>(?P<etb>\x17)|(?P<etx>\x03))$'

_MESSAGE_MATCHER = re.compile(_REGEXP, re.M | re.S)


def tokens( message ):
    m = _MESSAGE_MATCHER.match( message )
    if not m:
        raise MessageParseError(message)
    l = []
    for g in 'mode address tech_ack label ubi_dbi text end'.split():
        l += ( (g, m.group(g.upper())), )
    return dict(l)


def direction( ubi_dbi, decomposed_label ):
    direction = decomposed_label[0]
    if ubi_dbi is not None:
        if ubi_dbi.isdigit():
            direction = 'd'
        elif ubi_dbi.isalpha() or ubi_dbi == '\0':
            direction = 'u'
    return direction


_DTM_RE = r'^(?P<ORIGINATOR>.)(?P<MSN>..)(?P<BSC>.)' + \
          r'(?P<ALI>..)(?P<FN>....)(?P<FREE_TEXT>.*)$'
_DTM = re.compile(_DTM_RE, re.M | re.S)
_DL_MSG_KEYS = 'ORIGINATOR MSN BSC ALI FN FREE_TEXT'.split()

def dl_text(direction, text):
    ret = {}.fromkeys( (k.lower() for k in _DL_MSG_KEYS), 'UNKNOWN')
    ret['free_text'] = text
    if direction in 'db ':
        m = _DTM.match(text)
        if m is not None:
            ret.update( ((key.lower(), m.group(key)) for key in _DL_MSG_KEYS) )
    return ret


_WANTED_TOKENS = set("label ubi_dbi mode address tech_ack ali fn text end".split())
def parsed_message(raw_message):
    items = {}
    # Decompose raw_message.
    _tokens = tokens( raw_message )
    wanted_only = filter( lambda keyval: (keyval[0] in _WANTED_TOKENS), _tokens.iteritems() )
    items.update( wanted_only )

    # Remove the leading dots.
    items['address'] = items['address'].lstrip('.')

    # Label needs some special handling because it contains information about 
    # direction which is also contained in ubi_dbi, but ubi_dbi doesn't have to 
    # be present in the message according to old standard.
    label = items['label']
    ubi_dbi = items['ubi_dbi']
    _decomposed_label = decomposed_label( label )
    _direction = items['direction'] = direction( ubi_dbi, _decomposed_label )
    items['t_label'] = _decomposed_label[1]

    # Downlink messages has some little more of structured info - get it if 
    # possible!
    items.update( dl_text(_direction, items['text']) )
    return items


def t_mode( mode ):
    if '2' == mode:
        t_mode = "VHF Category A or Satellite Category A"
    elif 0x40 <= ord(mode) and ord(mode) <= 0x5D:
        t_mode = "VHF Category B, Ground station %s" % mode
    else:
        t_mode = "Unknown mode >%s<" % ord(mode)
    return t_mode


def t_free_text( free_text ):
    if None is free_text:
        free_text = " ||||||||||||| MESSAGE WITH NO TEXT ||||||||||||| " 
    return free_text


def t_tech_ack( tech_ack ):
    if tech_ack is None:
        tech_ack = ''
    if '\x15' == tech_ack: # Negative acknowledge.
        ret = "NAK"
    elif re.match(r"^[0-9a-zA-Z]$", tech_ack):
        ret = str(tech_ack)
    else:
        ret = "Unknown technical acknowledgment <%s>." \
            % tech_ack
    return ret


def t_originator( originator ):
    if originator is None:
        originator = ''
    default = 'Unknown originator <%s>' % originator
    return _ORIGINATOR_DICT.get(originator, default)
    

def t_direction( direction ):
    return {'d':'DWN', 'u':'UP', 'b': 'BTH' }.get(direction, 'UNK')


def t_msn( msn ):
    if msn is None:
        msn = ''
    return msn


def t_bsc( bsc ): 
    if bsc is None:
        bsc = ''
    return bsc

def t_end( end ): 
    if end is None:
        end = ''
    elif end == '\x17':
        end = "ETB"
    elif end == '\x03':
        end = "ETX"
    return end


class message(dict):
    """ An object representing the message. """

    def __init__(self, raw_text, utctime): # FIXME the name
        assert(isinstance(raw_text, basestring))
        self.original = raw_text
        self.update(parsed_message(raw_text))
        self.utcdatetime = utctime


    def __repr__(self):
        return self.original


    def __str__(self):
        m = []
        m += ( "Label:        %s" %               self['t_label'],     )
        m += ( "Flight:       %2s %4s" %        ( self['ali'],              
                                                  self['fn'],),   )
        m += ( "Address:      %s" %               self['address'],     )
        m += ( "Mode:         %s" %       t_mode( self['mode'] ),      )
        m += ( "Originator:   %s" % t_originator( self['originator']), )
        m += ( "Direction:    %s" %  t_direction( self['direction'] ), )
        m += ( "UDBI/TecAck:            %1s/%s" %
                     (     self['ubi_dbi'], 
               t_tech_ack( self['tech_ack'] ) ), )
        m += ( "MesSeqN/BlockSeqN/Suf:  %1s/%1s/%3s" \
                % ( t_msn( self['msn'] ), 
                    t_bsc( self['bsc'] ), 
                    t_end( self['end'] ) ), )
        m += ( '', )
        m += ( "Free text:", )
        m += ( t_free_text( self['free_text']), )
        return "\n".join(m)



class abstract_output_formatter:
    """ A base class for formatting the messages for the output. """

    def __init__(self, outfile):
        self.outfile = outfile

    def write_header(self):
        pass

    def write(self, m):
        raise NotImplementedError()
    
    def write_footer(self):
        pass



class csv_output_formater( abstract_output_formatter ):
    FIELD_NAMES=    'utc_year', 'utc_month', 'utc_day', \
                    'utc_hour', 'utc_minute', 'utc_second', \
                    'label', 't_label', 'ali', 'fn', 'address',  \
                    'mode', 't_mode', 'originator', 't_originator', \
                    'direction', 'ubi_dbi', 'msn', 't_msn', 'bsc', 't_bsc', \
                    'tech_ack', 'free_text', 't_end', 'raw'


    def __init__( self, outfile, writeheader=True):
        abstract_output_formatter.__init__( self, outfile )
        self.writer = DictWriter( outfile, self.FIELD_NAMES, 
                                    extrasaction='ignore' )
        self.writeheader = writeheader


    def write_header(self):
        if self.writeheader:
            # TODO With newer python writeheader can be used instead.
            print >>self.outfile, ",".join( self.FIELD_NAMES )


    def write(self, m):
        m = copy(m) # We may not want to the given message 
                    # object to get altered.
        utcdt = m.utcdatetime
        m['utc_year']   = str(utcdt.year)
        m['utc_month']  = str(utcdt.month)
        m['utc_day']    = str(utcdt.day)
        m['utc_hour']   = str(utcdt.hour)
        m['utc_minute'] = str(utcdt.minute)
        m['utc_second'] = str(utcdt.second+utcdt.microsecond * 1e-6)
    
        m['t_mode']         = str(t_mode(m['mode']))
        m['t_originator']   = str(t_originator(m['originator']))
        m['t_msn']          = str(t_msn(m['msn']))
        m['t_bsc']          = str(t_bsc(m['bsc']))
        m['t_end']          = str(t_end(m['end']))
        m['raw']            = m.original
        self.writer.writerow(m)



class hr_output_formater(abstract_output_formatter):
    def __init__(self, file):
        abstract_output_formatter.__init__(self, file )

    
    def write_header(self):
        print 


    def write(self, m):
        line = []
        line += ( "================================", )
        line += ( )
        line += ( "ACARS Message decoded at " \
                    + m.utcdatetime.strftime("[ %c UTC ]:"), )
        line += ( str(m), )
        line += ( )
        print >>self.outfile, "\n".join(line)


    def write_footer(self):
        print >>self.outfile, "=========== FINISHED ==========="



class raw_output_formater(abstract_output_formatter):
    def __init__(self, file, add_date=True):
        abstract_output_formatter.__init__(self, file )
        self.add_date = add_date
    

    def write(self, m):
        line = []
        if self.add_date:
            line += ( str( m.utcdatetime ) + "Z", )
        line += ( m.original, )
        self.outfile.write("\t".join(line) + "\n")
    
