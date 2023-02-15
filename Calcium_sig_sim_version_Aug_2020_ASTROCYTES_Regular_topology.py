# -*- coding: cp1252 -*-
#############################################################################################################################################
#############################################################################################################################################
##
##Calcium-Signalling-based Molecular Communications Simulator - CALCOMSIM (2023)
##
##Authors: Michael Taynnan Barros
##
##
###########################################################################################################################################
###########################################################################################################################################

from numpy import fabs
from random import uniform
from numpy.random import multinomial, exponential,random
from numpy import arange, array, empty,zeros,log
import sys
import time

##################################################################################
##Class responsible for the Gillespie method implementation
##################################################################################

class Model:
    def __init__(self,vnames,rates,inits, tmat,propensity,debug):
        '''
         * vnames: list of strings
         * rates: list of fixed rate parameters
         * inits: list of initial values of variables
         * propensity: list of lambda functions of the form:
            lambda r,ini: some function of rates ans inits.
        '''
        self.vn = vnames
        self.rates = rates
        self.inits = inits
        self.tm = tmat
        self.pv = propensity#[compile(eq,'errmsg','eval') for eq in propensity]
        self.pvl = len(self.pv) #length of propensity vector
        self.nvars = len(self.inits) #number of variables
        self.time=None
        self.series=None
        self.steps=0
        self.debug=debug

    def getStats(self):
        return self.time,self.series,self.steps

    def run(self,method='SSA', tmax=10, reps=1):
        self.res = zeros((tmax,self.nvars,reps),dtype=float)
        tvec = arange(tmax, dtype=int)
        if method =='SSA':
			rout2 = self.roulette
			for i in xrange(reps): ## I put the implementation of the GDA here in other to same time. Kept the function though. 21/03/2014
				'''
				Gillespie Direct algorithm
				'''
				ini = self.inits
				r = self.rates
				pvi = self.pv
				l=self.pvl
				pv = zeros(l,dtype=float)
				tm = self.tm

				tc = 0
				steps = 0
				self.res[0,:,i]= ini
				a0=1

				for i in xrange(l):
					pv[i] = pvi[i](ini,i)

				a0 = pv.sum() #sum of all transition probabilities
				tau = (-1/a0)*log(random())
				p='n'
				while p == 'n':
					###############################
					#####roulette function
					##############################
					limite = a0*random()
					aux = 0
					k = 0
					while((k<len(pv)) and (aux<limite)):
								aux = aux + pv[k]
								k = k + 1
					k = k - 1
					event = k
					p2 = uniform(0.0,1.0)
					sum1=0
					sum2=0
					for i in xrange(event):
						if l == i:
							sum1=sum2
						sum2+=pv[i]/a0
					if sum1 < p2 and sum2 >= p2:
						p='s'
				tc += tau
				if self.debug == 's':
					print 'ini :',ini
					print 'pv: ',pv
					print 'R :',event
					print 'Value: ', pv[event]
					print 'Release: ', pv[400]
					print 'RXACT: ', pv[416]
					print 'RXINACT', pv[417]
				pv2=pv

        elif method == 'SSAct':
            pass
        self.time=tvec
        self.series=self.res

        return tau,event, pv2

    def GSSA(self, tmax=50, round=0):
        '''
        Gillespie Direct algorithm
        '''
        ini = self.inits
        r = self.rates
        pvi = self.pv
        l=self.pvl
        pv = zeros(l,dtype=float)
        tm = self.tm

        tc = 0
        steps = 0
        self.res[0,:,round]= ini
        a0=1

        for i in xrange(l):
            pv[i] = pvi[i](ini,i)
        a0 = pv.sum() #sum of all transition probabilities
        tau = (-1/a0)*log(random())
        p='n'
        while p == 'n':
            event = self.roulette(pv,a0)
            p2 = uniform(0.0,1.0)
            sum1=0
            sum2=0
            for i in xrange(event):
                if l == i:
                    sum1=sum2
                sum2+=pv[i]/a0
            if sum1 < p2 and sum2 >= p2:
                p='s'
        tc += tau
        if self.debug == 's':
            print 'ini :',ini
            print 'pv: ',pv
            print 'R :',event
            print 'Value: ', pv[event]
            print 'Release: ', pv[400]
            print 'RXACT: ', pv[416]
            print 'RXINACT', pv[417]
        return tau, event, pv

    def roulette(self, popsel, popnumber):
        limite = popnumber*random()
        aux = 0
        k = 0
        while((k<len(popsel)) and (aux<limite)):
                    aux = aux + popsel[k]
                    k = k + 1
        k = k - 1
        return k

##################################################################################
#Main body of the simulator
##################################################################################

def main(destination, conc, n, m, x, txconc, deltaamp, tEnd):
    print '###############################################################'
    print 'Simulation of: 1) Distance from the Tx:', destination,'# of cells'
    print '2) Transmitter Concentration:', txconc,'nM' 
    print '3) Discrete 3D Tissue Size:', n, m
    print '4) Length of the cell:', x/10,'uM'
    print '###############################################################'

    ##################################################################################
    ##Variable Initialization
    ##################################################################################

    #### Probability of gap junction values astrocytes

    phl = [0.333333333333,0.0951756626745,0.0271812917035,0.00776661124715,0.00222288659849,0.000639921679991,0.000187410416804,5.81750667195e-05
    ,2.17596143238e-05,1.1436923158e-05,7.88454209682e-06,7.43738619183e-06,7.37970057786e-06,7.29603316347e-06,7.27478942971e-06,7.26006289992e-06
    ,7.26084787208e-06,7.26080132601e-06,7.26061054996e-06,7.26081361742e-06,7.26079620991e-06,7.26072365567e-06,7.26058079345e-06,7.26074725419e-06
    ,7.26087576894e-06,7.26073008288e-06,7.26061194028e-06,7.26074336727e-06,7.26101528686e-06,7.26081974085e-06,7.26091667847e-06,7.26058059059e-06
    ,7.26084014577e-06,7.2610063969e-06,7.26069065682e-06,7.26083092741e-06,7.26076153595e-06,7.26071756287e-06,7.26092535023e-06,7.26076421324e-06
    ,7.26060026219e-06,7.26075209967e-06,7.26093367537e-06,7.26073986493e-06,7.26039032094e-06,7.26091299989e-06,7.26077756319e-06,7.26071491915e-06
    ,7.2607710224e-06,7.26082337127e-06]

    plh = [0.333333333333,0.706523083189,0.825908352326,0.86117728508,0.871356970354,0.874272895353,0.875107602476,0.875346523333,0.875414979813
    ,0.875433066086,0.875439533589,0.875437420808,0.875439605921,0.875443525702,0.875437685437,0.875440602592,0.875441218644,0.875440938251
    ,0.875438148928,0.875441400815,0.875441147657,0.875440025483,0.875437801975,0.875440355641,0.875442338151,0.875440081279,0.87543825197
    ,0.875440285023,0.87544449347,0.875441466845,0.875442967173,0.875437765288,0.875441782611,0.875444355802,0.875439468868,0.875441639938
    ,0.875440565916,0.875439885313,0.875443101387,0.875440607355,0.875438069767,0.875440419864,0.875443230241,0.875440230498,0.875434820356
    ,0.875442910232,0.875440813981,0.875439844395,0.875440712745,0.875441522986]

    phh = [0.333333333333,0.198301254137,0.14691035597,0.131056103673,0.126420143048,0.125087182967,0.124704987107,0.1245953016,0.124563260573
    ,0.124555496991,0.124552581869,0.124555141806,0.124553014379,0.124549178265,0.124555039774,0.124552137345,0.124551520508,0.124551800947
    ,0.124554590462,0.124551338371,0.124551591547,0.124552713793,0.124554937445,0.124552383612,0.124550400973,0.124552657991,0.124554487418
    ,0.124552454233,0.124548245514,0.124551272335,0.12454977191,0.124554974131,0.124550956549,0.124548383192,0.124553270441,0.124551099232
    ,0.124552173322,0.124552853969,0.124549637687,0.124552131881,0.124554669633,0.124552319384,0.124549508825,0.124552508762,0.124557919254
    ,0.124549828855,0.124551925241,0.12455289489,0.124552026484,0.124551216191]

    D = 122500
    l = 0.5
    space = 10
    laticeset = [0] * int(space/l) * n * m
    t = 0
    #tEnd = 4  # Total runtime
    tini = 0
    tend = 2
    C = []
    C2 = []
    T = []
    numberofrec = 29
    numberofvar = 30
    ini = [100,0,0.05,15,0.1,2.02,0.1,4000,50,1.5,50,0.3,8,0.05,0.15,0.15,0.1,0.1,4,0,0,0.5,D,(3.1416*((float(x)/2.0)**2)),0.0006,2.5,2.2,phh[0],phl[0],plh[0]] * int(space/l) * n * m
    vsars = ['0v1','1Y','2vin','3VM2','4C','5n','6K2','7VM3','8kout','9S','10kf','11kp','12kdeg','13vp','14kcaaa','15kcai','16kip3','17Z'
    ,'18q','19W','20A','21kia','22D','23l','24K','25ka','26m','27phh','28phl','29plh'] * int(space/l) * n * m
    prop = []
    freq = [0]*(numberofrec)
    freq2 = [0]*int(space/l) * n * m
    ALPHA = 0.01
    q = 4
    f=0
    states=0
    stater=1
    liststates=[]
    liststater=[]
    ini[(len(laticeset)/2)*numberofvar+1]=txconc
    ini[(len(laticeset)/2)*numberofvar+4]=deltaamp
    dest = destination
    ini[(len(laticeset)/2+dest)*numberofvar+19]=conc
    destinationid = (len(laticeset)/2+dest)*numberofvar+4
    listdiffTx = [0]*4
    listdiffRx = [0]*4
    debug='n'
    p=[]
    p2=[]
    p4=[]
    p5=[]
    lpn=[]
    maxGainSamples = 10
    GainTrans = []
    GainRec = []

    loslist = []
    loslistappend = loslist.append

    for c in xrange((len(laticeset)/2)+1,(len(laticeset)/2)+dest):
        loslistappend(c)

    a = 1.05 # linear coeficient

    props1=0
    props0=0
    propr1=0
    propr0=0

    probr1s0=0
    probr0s0=0
    probr1s1=0
    probr0s1=0

    propappend = prop.append
    propextend = prop.extend

    ##################################################################################
    ## Initialization of reaction functions
    ##################################################################################

    for tant in xrange(0,len(ini),numberofvar):
        # reactions
        propextend([lambda ini,i:ini[(i/numberofrec)*numberofvar+0]*ini[(i/numberofrec)*numberofvar+1]
        ,lambda ini,i:ini[(i/numberofrec)*numberofvar+2]
        ,lambda ini,i:4*ini[(i/numberofrec)*numberofvar+7]*((ini[(i/numberofrec)*numberofvar+14]**ini[(i/numberofrec)*numberofvar+5]*ini[(i/numberofrec)*numberofvar+4]**ini[(i/numberofrec)*numberofvar+5])/((ini[(i/numberofrec)*numberofvar+4]**ini[(i/numberofrec)*numberofvar+5]+ini[(i/numberofrec)*numberofvar+14]**ini[(i/numberofrec)*numberofvar+5])*(ini[(i/numberofrec)*numberofvar+4]**ini[(i/numberofrec)*numberofvar+5]+ini[(i/numberofrec)*numberofvar+15]**ini[(i/numberofrec)*numberofvar+5])))*
					((ini[(i/numberofrec)*numberofvar+17]**ini[(i/numberofrec)*numberofvar+26])/(ini[(i/numberofrec)*numberofvar+17]**ini[(i/numberofrec)*numberofvar+26] + ini[(i/numberofrec)*numberofvar+16]**ini[(i/numberofrec)*numberofvar+26]))*(ini[(i/numberofrec)*numberofvar+9]-ini[(i/numberofrec)*numberofvar+4])
        ,lambda ini,i:ini[(i/numberofrec)*numberofvar+3]*((ini[(i/numberofrec)*numberofvar+4]**2)/(ini[(i/numberofrec)*numberofvar+4]**2+ini[(i/numberofrec)*numberofvar+6]**2))
        ,lambda ini,i:ini[(i/numberofrec)*numberofvar+10]*ini[(i/numberofrec)*numberofvar+9]
        ,lambda ini,i:ini[(i/numberofrec)*numberofvar+10]*ini[(i/numberofrec)*numberofvar+4]
        ,lambda ini,i:ini[(i/numberofrec)*numberofvar+8]*ini[(i/numberofrec)*numberofvar+4]
        ,lambda ini,i:ini[(i/numberofrec)*numberofvar+13]*((ini[(i/numberofrec)*numberofvar+2]**2)/(ini[(i/numberofrec)*numberofvar+4]**2+ini[(i/numberofrec)*numberofvar+11]**2))
        ,lambda ini,i:ini[(i/numberofrec)*numberofvar+12]*ini[(i/numberofrec)*numberofvar+17]])

        ## receiver reactions
        if tant == (len(laticeset)/2+dest)*numberofvar:
            print 'Communication System is ready'
            propappend(lambda ini,i:ini[(i/numberofrec)*numberofvar+25]*((ini[(i/numberofrec)*numberofvar+4]**ini[(i/numberofrec)*numberofvar+18])/(ini[(i/numberofrec)*numberofvar+24]**ini[(i/numberofrec)*numberofvar+18]+ini[(i/numberofrec)*numberofvar+4]**ini[(i/numberofrec)*numberofvar+18]))
                    *(ini[(i/numberofrec)*numberofvar+19]-conc*ini[(i/numberofrec)*numberofvar+20]))
            propappend(lambda ini,i:ini[(i/numberofrec)*numberofvar+21]*conc*ini[(i/numberofrec)*numberofvar+20])
        else:
            propappend(lambda ini,i:ini[(i/numberofrec)*numberofvar+25]*((ini[(i/numberofrec)*numberofvar+4]**ini[(i/numberofrec)*numberofvar+18])/(ini[(i/numberofrec)*numberofvar+24]**ini[(i/numberofrec)*numberofvar+18]+ini[(i/numberofrec)*numberofvar+4]**ini[(i/numberofrec)*numberofvar+18]))
                    *(ini[(i/numberofrec)*numberofvar+19]-0*ini[(i/numberofrec)*numberofvar+20]))
            propappend(lambda ini,i:ini[(i/numberofrec)*numberofvar+21]*0*ini[(i/numberofrec)*numberofvar+20])

        ## DIFFUSION FOR REACTIONS

        ## C_x,y,z -> C_x+1,y,z

        ##phh
        if tant+numberofvar < len(ini):
            propappend(lambda ini,i:(((ini[(i/numberofrec)*numberofvar+22]/(ini[(i/numberofrec)*numberofvar+23]))*fabs(ini[(i/numberofrec)*numberofvar+4]-ini[(i/numberofrec)*numberofvar+numberofvar+4]))*ini[(i/numberofrec)*numberofvar+27]) if ini[(i/numberofrec)*numberofvar+4] > ini[(i/numberofrec)*numberofvar+numberofvar+4] else 0)
        else:
            propappend(lambda ini,i:0)

        ##phl
        if tant+numberofvar < len(ini):
            propappend(lambda ini,i:(((ini[(i/numberofrec)*numberofvar+22]/(ini[(i/numberofrec)*numberofvar+23]))*fabs(ini[(i/numberofrec)*numberofvar+4]-ini[(i/numberofrec)*numberofvar+numberofvar+4]))*ini[(i/numberofrec)*numberofvar+28]) if ini[(i/numberofrec)*numberofvar+4] > ini[(i/numberofrec)*numberofvar+numberofvar+4] else 0)
        else:
            propappend(lambda ini,i:0)

        ##plh
        if tant+numberofvar < len(ini):
            propappend(lambda ini,i:(((ini[(i/numberofrec)*numberofvar+22]/(ini[(i/numberofrec)*numberofvar+23]))*fabs(ini[(i/numberofrec)*numberofvar+4]-ini[(i/numberofrec)*numberofvar+numberofvar+4]))*ini[(i/numberofrec)*numberofvar+29]) if ini[(i/numberofrec)*numberofvar+4] > ini[(i/numberofrec)*numberofvar+numberofvar+4] else 0)
        else:
            propappend(lambda ini,i:0)


        ## C_x,y,z -> C_x-1,y,z

        ##phh
        if tant-numberofvar > 0:
            propappend(lambda ini,i:(((ini[(i/numberofrec)*numberofvar+22]/(ini[(i/numberofrec)*numberofvar+23]))*fabs(ini[(i/numberofrec)*numberofvar+4]-ini[(i/numberofrec)*numberofvar-numberofvar+4]))*ini[(i/numberofrec)*numberofvar+27]) if ini[(i/numberofrec)*numberofvar+4] > ini[(i/numberofrec)*numberofvar-numberofvar+4] else 0)
        else:
            propappend(lambda ini,i:0)

        ##phl
        if tant-numberofvar > 0:
            propappend(lambda ini,i:(((ini[(i/numberofrec)*numberofvar+22]/(ini[(i/numberofrec)*numberofvar+23]))*fabs(ini[(i/numberofrec)*numberofvar+4]-ini[(i/numberofrec)*numberofvar-numberofvar+4]))*ini[(i/numberofrec)*numberofvar+28]) if ini[(i/numberofrec)*numberofvar+4] > ini[(i/numberofrec)*numberofvar-numberofvar+4] else 0)
        else:
            propappend(lambda ini,i:0)

        ##plh
        if tant-numberofvar > 0:
            propappend(lambda ini,i:(((ini[(i/numberofrec)*numberofvar+22]/(ini[(i/numberofrec)*numberofvar+23]))*fabs(ini[(i/numberofrec)*numberofvar+4]-ini[(i/numberofrec)*numberofvar-numberofvar+4]))*ini[(i/numberofrec)*numberofvar+29]) if ini[(i/numberofrec)*numberofvar+4] > ini[(i/numberofrec)*numberofvar-numberofvar+4] else 0)
        else:
            propappend(lambda ini,i:0)

        ## C_x,y,z -> C_x,y+1,z

        ##phh
        if tant+(numberofvar * int(space/l)) < len(ini):
            propappend(lambda ini,i:(((ini[(i/numberofrec)*numberofvar+22]/(ini[(i/numberofrec)*numberofvar+23]))*fabs(ini[(i/numberofrec)*numberofvar+4]-ini[(i/numberofrec+int(space/l))*numberofvar+4]))*ini[(i/numberofrec)*numberofvar+27]) if ini[(i/numberofrec)*numberofvar+4] > ini[(i/numberofrec+int(space/l))*numberofvar+4] else 0)
        else:
            propappend(lambda ini,i:0)

        ##phl
        if tant+(numberofvar * int(space/l)) < len(ini):
            propappend(lambda ini,i:(((ini[(i/numberofrec)*numberofvar+22]/(ini[(i/numberofrec)*numberofvar+23]))*fabs(ini[(i/numberofrec)*numberofvar+4]-ini[(i/numberofrec+int(space/l))*numberofvar+4]))*ini[(i/numberofrec)*numberofvar+28]) if ini[(i/numberofrec)*numberofvar+4] > ini[(i/numberofrec+int(space/l))*numberofvar+4] else 0)
        else:
            propappend(lambda ini,i:0)

        #plh
        if tant+(numberofvar * int(space/l)) < len(ini):
            propappend(lambda ini,i:(((ini[(i/numberofrec)*numberofvar+22]/(ini[(i/numberofrec)*numberofvar+23]))*fabs(ini[(i/numberofrec)*numberofvar+4]-ini[(i/numberofrec+int(space/l))*numberofvar+4]))*ini[(i/numberofrec)*numberofvar+29]) if ini[(i/numberofrec)*numberofvar+4] > ini[(i/numberofrec+int(space/l))*numberofvar+4] else 0)
        else:
            propappend(lambda ini,i:0)

        ## C_x,y,z -> C_x,y-1,z

        ##phh
        if tant-(numberofvar * int(space/l)) > 0:
            propappend(lambda ini,i:(((ini[(i/numberofrec)*numberofvar+22]/(ini[(i/numberofrec)*numberofvar+23]))*fabs(ini[(i/numberofrec)*numberofvar+4]-ini[(i/numberofrec-int(space/l))*numberofvar+4]))*ini[(i/numberofrec)*numberofvar+27]) if ini[(i/numberofrec)*numberofvar+4] > ini[(i/numberofrec-int(space/l))*numberofvar+4] else 0)
        else:
            propappend(lambda ini,i:0)

        #phl
        if tant-(numberofvar * int(space/l)) > 0:
            propappend(lambda ini,i:(((ini[(i/numberofrec)*numberofvar+22]/(ini[(i/numberofrec)*numberofvar+23]))*fabs(ini[(i/numberofrec)*numberofvar+4]-ini[(i/numberofrec-int(space/l))*numberofvar+4]))*ini[(i/numberofrec)*numberofvar+28]) if ini[(i/numberofrec)*numberofvar+4] > ini[(i/numberofrec-int(space/l))*numberofvar+4] else 0)
        else:
            propappend(lambda ini,i:0)

        #plh
        if tant-(numberofvar * int(space/l)) > 0:
            propappend(lambda ini,i:(((ini[(i/numberofrec)*numberofvar+22]/(ini[(i/numberofrec)*numberofvar+23]))*fabs(ini[(i/numberofrec)*numberofvar+4]-ini[(i/numberofrec-int(space/l))*numberofvar+4]))*ini[(i/numberofrec)*numberofvar+29]) if ini[(i/numberofrec)*numberofvar+4] > ini[(i/numberofrec-int(space/l))*numberofvar+4] else 0)
        else:
            propappend(lambda ini,i:0)

        ## C_x,y,z -> C_x,y,z+1

        #phh
        if tant+(numberofvar * int(space/l) * n) < len(ini):
            propappend(lambda ini,i:(((ini[(i/numberofrec)*numberofvar+22]/(ini[(i/numberofrec)*numberofvar+23]))*fabs(ini[(i/numberofrec)*numberofvar+4]-ini[(i/numberofrec+(int(space/l)*n))*numberofvar+4]))*ini[(i/numberofrec)*numberofvar+27]) if ini[(i/numberofrec)*numberofvar+4] > ini[(i/numberofrec+(int(space/l)*n))*numberofvar+4] else 0)
        else:
            propappend(lambda ini,i:0)

        #phl
        if tant+(numberofvar * int(space/l) * n) < len(ini):
            propappend(lambda ini,i:(((ini[(i/numberofrec)*numberofvar+22]/(ini[(i/numberofrec)*numberofvar+23]))*fabs(ini[(i/numberofrec)*numberofvar+4]-ini[(i/numberofrec+(int(space/l)*n))*numberofvar+4]))*ini[(i/numberofrec)*numberofvar+28]) if ini[(i/numberofrec)*numberofvar+4] > ini[(i/numberofrec+(int(space/l)*n))*numberofvar+4] else 0)
        else:
            propappend(lambda ini,i:0)

        #plh
        if tant+(numberofvar * int(space/l) * n) < len(ini):
            propappend(lambda ini,i:(((ini[(i/numberofrec)*numberofvar+22]/(ini[(i/numberofrec)*numberofvar+23]))*fabs(ini[(i/numberofrec)*numberofvar+4]-ini[(i/numberofrec+(int(space/l)*n))*numberofvar+4]))*ini[(i/numberofrec)*numberofvar+29]) if ini[(i/numberofrec)*numberofvar+4] > ini[(i/numberofrec+(int(space/l)*n))*numberofvar+4] else 0)
        else:
            propappend(lambda ini,i:0)

        ## C_x,y,z -> C_x,y,z-1

        #phh
        if tant-(numberofvar * int(space/l) * n) > 0:
            propappend(lambda ini,i:(((ini[(i/numberofrec)*numberofvar+22]/(ini[(i/numberofrec)*numberofvar+23]))*fabs(ini[(i/numberofrec)*numberofvar+4]-ini[(i/numberofrec-(int(space/l)*n))*numberofvar+4]))*ini[(i/numberofrec)*numberofvar+27]) if ini[(i/numberofrec)*numberofvar+4] > ini[(i/numberofrec-(int(space/l)*n))*numberofvar+4] else 0)
        else:
            propappend(lambda ini,i:0)

        #phl
        if tant-(numberofvar * int(space/l) * n) > 0:
            propappend(lambda ini,i:(((ini[(i/numberofrec)*numberofvar+22]/(ini[(i/numberofrec)*numberofvar+23]))*fabs(ini[(i/numberofrec)*numberofvar+4]-ini[(i/numberofrec-(int(space/l)*n))*numberofvar+4]))*ini[(i/numberofrec)*numberofvar+28]) if ini[(i/numberofrec)*numberofvar+4] > ini[(i/numberofrec-(int(space/l)*n))*numberofvar+4] else 0)
        else:
            propappend(lambda ini,i:0)

        #plh
        if tant-(numberofvar * int(space/l) * n) > 0:
            propappend(lambda ini,i:(((ini[(i/numberofrec)*numberofvar+22]/(ini[(i/numberofrec)*numberofvar+23]))*fabs(ini[(i/numberofrec)*numberofvar+4]-ini[(i/numberofrec-(int(space/l)*n))*numberofvar+4]))*ini[(i/numberofrec)*numberofvar+29]) if ini[(i/numberofrec)*numberofvar+4] > ini[(i/numberofrec-(int(space/l)*n))*numberofvar+4] else 0)
        else:
            propappend(lambda ini,i:0)

    flag = True

    ##################################################################################
    ## Initialization of the Simulator
    ##################################################################################
    print 'Simulation has been started. Wait for the results\n'

    # setup toolbar
    progress = 0
    toolbar_width = 50
    sys.stdout.write("PROGRESS: [%s]" % (" " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['
    sys.stdout.flush()

    while t <= tEnd:

        r = [.001,.1]
        tm = 2

        M = Model(vnames = vars,rates = r,inits=ini, tmat=tm,propensity=prop,debug=debug)
        tau,event,pv = M.run(tmax=200,reps=1)
        t = t + tau*1000
        reaction = event
        i=reaction/numberofrec
        x=(reaction-(reaction/numberofrec)*numberofrec)+1

        # This flag is necessary to perform the first iteration
        if (flag):
            flag = False
            t_ant = t

        if int(t) != f:
            f = int(t)
            if (tEnd < toolbar_width):
                progress += toolbar_width / tEnd
                sys.stdout.write("%s" % ("#" * (toolbar_width / tEnd)))
                sys.stdout.flush()
                if (f == tEnd and progress != toolbar_width):
                    sys.stdout.write("%s" % ("#" * (toolbar_width - progress)))
                    sys.stdout.flush()
            else:
                if (f % (tEnd / toolbar_width) == 0 and f != 0):
                    # update the bar
                    sys.stdout.write("#")
                    sys.stdout.flush()
                    if (f == tEnd and progress != toolbar_width):
                        sys.stdout.write("%s" % ("#" * (toolbar_width - progress)))
                        sys.stdout.flush()

            itme = 0

            ##update phh, phl and plh
            for tant in xrange(0,len(ini),numberofvar):
				if int(t) < len(phh):
					ini[tant+27] = phh[int(t)]
					ini[tant+28] = phl[int(t)]
					ini[tant+29] = plh[int(t)]
				else:
					ini[tant+27] = phh[len(phh)-1]
					ini[tant+28] = phl[len(phh)-1]
					ini[tant+29] = plh[len(phh)-1]

    ##################################################################################
    ## Reactions Execution
    ##################################################################################

        if int(t) == tini:
            states=1
        if int (t) == tend:
            states=0
        if x == 2:
            ini[i*numberofvar+4]+=ALPHA
        if x == 3:
            ini[i*numberofvar+4]+=ALPHA
            ini[i*numberofvar+9]-=ALPHA
            if ini[i*numberofvar+9] < 0 :
                ini[i*numberofvar+9] = ini[i*numberofvar+9]*(-1)
        if x == 4:
            ini[i*numberofvar+4]-=ALPHA
            if ini[i*numberofvar+4] < 0:
                ini[i*numberofvar+4] = ini[i*numberofvar+4]*(-1)
            ini[i*numberofvar+9]+=ALPHA
        if x == 5:
            ini[i*numberofvar+4]+=ALPHA
            ini[i*numberofvar+8]-=ALPHA
            if ini[i*numberofvar+9] < 0:
                ini[i*numberofvar+9] = ini[i*numberofvar+9]*(-1)
        if x == 6:
            ini[i*numberofvar+4]-=ALPHA
            ini[i*numberofvar+9]+=ALPHA
            if ini[i*numberofvar+4] < 0:
                ini[i*numberofvar+4] = ini[i*numberofvar+4]*(-1)

        if x == 7:
            ini[i*numberofvar+4]-=ALPHA
            if ini[i*numberofvar+4] < 0:
                ini[i*numberofvar+4] = ini[i*numberofvar+4]*(-1)

        if x == 8:
            ini[i*numberofvar+17]-=ALPHA

        if x == 9:
            ini[i*numberofvar+17]-=ALPHA
            if ini[i*numberofvar+17] < 0:
                ini[i*numberofvar+17] = ini[i*numberofvar+17]*(-1)

        if x == 10 :
            ini[i*numberofvar+4]-=q*ALPHA
            ini[i*numberofvar+19]+=q*ALPHA
            if i == len(laticeset)/2+dest:
                stater=1
                ini[(len(laticeset)/2+dest)*numberofvar+20]+=1
            if ini[i*numberofvar+4] < 0:
                ini[i*numberofvar+4] = ini[i*numberofvar+4]*(-1)

        if x == 11 :
            ini[i*numberofvar+4]+=q*ALPHA
            ini[i*numberofvar+19]-=q*ALPHA
            if i == len(laticeset)/2+dest:
                stater=0
                ini[(len(laticeset)/2+dest)*numberofvar+20]-=1
            if ini[i*numberofvar+19] < 0:
                ini[i*numberofvar+19] = ini[i*numberofvar+19]*(-1)

        ## DIFFUSION REACTIONS START HERE

        if x == 12:
            ini[i*numberofvar+4]-=ALPHA
            if ini[i*numberofvar+4] < 0:
                ini[i*numberofvar+4] = ini[i*numberofvar+4]*(-1)
            if i < len(laticeset)-1:
                ini[(i+1)*numberofvar+4]+=ALPHA
            freq2[i]=freq2[i]+1
            if i == len(laticeset)/2:
                listdiffTx[0] = listdiffTx[0]+1
            if i == len(laticeset)/2+dest:
                listdiffRx[0] = listdiffRx[0]+1

        if x == 15:
            ini[i*numberofvar+4]-=ALPHA
            if ini[i*numberofvar+4] < 0:
                ini[i*numberofvar+4] = ini[i*numberofvar+4]*(-1)
            if i != 0:
                ini[(i-1)*numberofvar+4]+=ALPHA
            freq2[i]=freq2[i]+1
            if i == len(laticeset)/2:
                listdiffTx[1] = listdiffTx[1]+1
            if i == len(laticeset)/2+dest:
                listdiffRx[1] = listdiffRx[1]+1

        if x == 18:
            ini[i*numberofvar+4]-=ALPHA
            if ini[i*numberofvar+4] < 0:
                ini[i*numberofvar+4] = ini[i*numberofvar+4]*(-1)
            if i < len(laticeset)-1:
                ini[(i+int(space/l))*numberofvar+4]+=ALPHA
            freq2[i]=freq2[i]+1
            if i == len(laticeset)/2:
                listdiffTx[2] = listdiffTx[2]+1
            if i == len(laticeset)/2+dest:
                listdiffRx[2] = listdiffRx[2]+1

        if x == 21:
            ini[i*numberofvar+4]-=ALPHA
            if ini[i*numberofvar+4] < 0:
                ini[i*numberofvar+4] = ini[i*numberofvar+4]*(-1)
            if i != 0:
                ini[(i-int(space/l))*numberofvar+4]+=ALPHA
            freq2[i]=freq2[i]+1
            if i == len(laticeset)/2:
                listdiffTx[2] = listdiffTx[3]+1
            if i == len(laticeset)/2+dest:
                listdiffRx[2] = listdiffRx[3]+1

        if x == 24:
            ini[i*numberofvar+4]-=ALPHA
            if ini[i*numberofvar+4] < 0:
                ini[i*numberofvar+4] = ini[i*numberofvar+4]*(-1)
            if i < len(laticeset)-1:
                ini[(i+(int(space/l)*n))*numberofvar+4]+=ALPHA
            freq2[i]=freq2[i]+1
            if i == len(laticeset)/2:
                listdiffTx[2] = listdiffTx[2]+1
            if i == len(laticeset)/2+dest:
                listdiffRx[2] = listdiffRx[2]+1

        if x == 27:
            ini[i*numberofvar+4]-=ALPHA
            if ini[i*numberofvar+4] < 0:
                ini[i*numberofvar+4] = ini[i*numberofvar+4]*(-1)
            if i != 0:
                ini[(i-(int(space/l)*n))*numberofvar+4]+=ALPHA
            freq2[i]=freq2[i]+1
            if i == len(laticeset)/2:
                listdiffTx[2] = listdiffTx[3]+1
            if i == len(laticeset)/2+dest:
                listdiffRx[2] = listdiffRx[3]+1


        freq[x-1]= freq[x-1] + 1

        C.append(ini[(len(laticeset)/2+dest)*numberofvar+4])
        C2.append(ini[(len(laticeset)/2)*numberofvar+4])
        T.append(t)

        if t <= maxGainSamples:
			GainTrans.append(ini[(len(laticeset)/2)*numberofvar+4])
			GainRec.append(ini[(len(laticeset)/2+dest)*numberofvar+4])

        liststates.append(states)
        liststater.append(stater)

        if debug == 's':
            cont = raw_input('\nContinua:')
            if cont is not 's':
                break
            else:
                continue

        if states == 1 and stater == 1 and (t >= tini and t <= tend):
            p.append(1)
        else:
            p.append(0)

        if states == 0 and stater == 1 and (t >= tini and t <= tend):
            p2.append(1)
        else:
            p2.append(0)

        if float(liststater.count(1))/float(len(liststater))*(ini[(len(laticeset)/2+dest)*numberofvar+4]) != 0:
            p4.append( ( (ini[(len(laticeset)/2+dest)*numberofvar+4]) / (ini[(len(laticeset)/2)*numberofvar+1]*ini[(len(laticeset)/2)*numberofvar])) )
        else:
            p4.append(0)

        props1=float(sum(liststates))/float(len(liststates))
        props0=1-props1
        propr1=float(liststater.count(1))/float(len(liststater))
        propr0=1-propr1

        if propr0 != 0 and propr1 != 0:
            probr1s0=(props0*propr1)/propr1
            probr0s0=1-probr1s0
            probr1s1=(props1*propr1)/propr1
            probr0s1=1-probr1s1


        p5.append(propr1*probr1s1)

        if p4[-1] >= 0.2 :
                     lpn.append(1)
        else:
			lpn.append(0)


    print '\n\n###### Number of Diffusions Reactions #####'
    print 'Tx: ',listdiffTx
    print 'Rx: ',listdiffRx

    print ""
    print "###### Mutual Information #####"
    print ""

    px1 = float(sum(liststates))/float(len(liststates))
    px0 = 1 - px1

    py1 = float(sum(p4))/float(len(p4))
    py0 = 1 - py1

    if py1==0:
        py1=0.00000000001
        py0 =0.9999999999

    py1x0 = float(sum(p2))/float(len(p2))
    py0x0 = 1 - py1x0

    if py1x0==0:
        py1x0=0.00000000001
        py0x0 =0.9999999999
    py1x1 = py1+py1x0-(py1*py1x0)
    py0x1 = 1 - py1x1

    from math import log

    ibits= px0*py0x0*(log(py0x0/(py0))/log(2))+px1*py0x1*(log(py0x1/(py0))/log(2))+px0*py1x0*(log(py1x0/(py1))/log(2))+px1*py1x1*(log(py1x1/(py1))/log(2))

    print ibits 


##################################################################################
##Calling up the simulations
##################################################################################

main(1,100,3,3,5,2,0.5,4)
