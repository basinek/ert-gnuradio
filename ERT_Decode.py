#!/usr/bin/env python
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
## OOT GNURadio Module: Decoder for ERT standard. v1.0: August 30, 2015
#  A.G. Klein and N. Conroy and K. Basinet
#  Description:  -Decodes ERT gas meter data from .bin file and displays
#                 meter ID, meter type, physical tamper flag, encoder tamper flag,
#                 consumption value, and the elapsed time since the last loop.
#                 Code is based on similary ERT decoder for MATLAB written by 
#                 A.G. Klein and N. Conroy and K. Basinet.
#                -Inputs: Complex vector with length of data block size
#                -Outputs: Complex vector with length of data block size
#  Dependencies: -Requires Numpy, GNURadio and custom function polynomialDivision
#  ------------------------------------------------------------------------
import numpy as np
from polynomialDivision import polynomialDivision
import numpy
from gnuradio import gr

# Parameters and constants
JMP = 30                            # Number of samples to jump over each iteration 
DataRate = 16384                    # Data rate for determining symbol period
SMPRT = 2392064                     # RTL-SDR Sample Rate
BLOCKSIZE = 18688                   # RTL-SDR Samples per frame
SP = numpy.int16(SMPRT/DataRate)          # Nominal symbol period (in # samples)
BCH_POLY = [1,0,1,1,0,1,1,1,1,0,1,1,0,0,0,1,1] # BCH generator polynomial coefficients from ERT standard
PREAMBLE = [1,1,1,1,1,0,0,1,0,1,0,1,0,0,1,1,0,0,0,0,0]  #From ERT standard, includes sync bit

#"Macro" for converting binary number (each digit as number in a list) to decimal integer
def bin2dec(bin_list):
    bin_list = [int(b) for b in bin_list]
    return int(''.join(str(c) for c in bin_list),2)

class ERT_Decode(gr.sync_block):
    """
    docstring for block blk1
    """

    def __init__(self, multiple):
        self.multiple = multiple
        self.data_buffer = []
        gr.sync_block.__init__(self,
            name="ERT Decoder",
            in_sig=[(numpy.uint8,18688)],
            out_sig=[(numpy.uint8,18688)])

    def work(self, input_items, output_items):
        output_items[0][0] = input_items[0][0]
        out = output_items[0][0]
        if (len(self.data_buffer)>(BLOCKSIZE*30)):
            dat=self.data_buffer-127
            self.data_buffer = []
            s = dat[1:(len(dat)-1):2]+1j*dat[2:(len(dat)-1):2]
            
            #Preallocate buffer space
            zbuff = numpy.zeros(BLOCKSIZE)
            softbits = numpy.zeros(96)
            bits = numpy.zeros(96)
            cnt = 0 #Decoded message counter
            block_index = 0
            while block_index < len(s)-BLOCKSIZE+JMP:
                i=0 # Counter for sample feeding
                zbuff = s[block_index:block_index+(BLOCKSIZE-1)] #Grab 18688 samples from file, store them in buffer
                buff = numpy.int32(numpy.real(zbuff))**2+numpy.int32(numpy.imag(zbuff))**2 #Cheap absolute value of buffer
                while i < BLOCKSIZE-(96*SP):
                    cu = numpy.cumsum(buff[i:(i+96*SP)])
                    softbits = (2*cu[(SP/2)+1:(95*SP)+(SP/2)+1:SP])-cu[1:(95*SP)+1:SP]-cu[SP+1:(95*SP)+SP+1:SP];
                    for n in range(len(softbits)): #List with '1' where corresponding index in softbits is positive
                        if softbits[n] > 0:
                            bits[n] = 1
                        else:
                            bits[n] = 0
                    #Check if preamble is correct and parse data
                    if numpy.array_equal(bits[0:len(PREAMBLE)],PREAMBLE):
                        #BCH processing
                        dc = numpy.concatenate([numpy.zeros(180),bits[21:96]])
                        if 1==1: #polynomialDivision(BCH_POLY,bits[21:96])[0] == 0:
                            #BCH passed
                            i = i+(96*SP)-JMP
                            cnt = cnt+1
                            #Separate BCH decoded blocks
                            dc_id = numpy.concatenate([dc[180:182],dc[215:239]])
                            SCM_ID = numpy.concatenate([bits[21:23],bits[55:79]])
                            dc_phy_tmp = dc[183:185]
                            dc_ert_type = dc[185:189]
                            dc_enc_tmp = dc[189:191]
                            dc_consump = dc[191:215]
                            #Convert to decimal
                            dc_id = bin2dec(dc_id)
                            dc_phy_tmp = bin2dec(dc_phy_tmp)
                            dc_ert_type = bin2dec(dc_ert_type)
                            dc_enc_tmp = bin2dec(dc_enc_tmp)
                            dc_consump = bin2dec(dc_consump)
                            #Print decoded output
                            print("Decoded Meter ID: %u" %dc_id)
                            print("Decoded Meter Type: %u" %dc_ert_type)
                            print("Decoded Physical Tamper: %u" %dc_phy_tmp)
                            print("Decoded Encoder Tamper: %u" %dc_enc_tmp)
                            print("Decoded Consumption: %u \n" %dc_consump)
                    i = i+JMP
                block_index = block_index+(JMP*96)
        else:
            self.data_buffer = numpy.concatenate([self.data_buffer,output_items[0][0]])
        return len(output_items[0])
#Done
