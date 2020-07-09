# -*- coding: cp1252 -*-
#############################################################################################################################################
#############################################################################################################################################
##
##Calcium signaling simulator (Version April 2014)
##
##Author: Michael Taynnan Barros (michael.taob@gmail.com)
##
##				New features: ASTROCYTES.
##
##
##This software is part of the author's Ph.D. and it is exclusively made for research purposes.
##The use of such software should be authorized by its author.
##
##We used Guillespie exact stochastic method alongside with the mininal model for calcium wave oscillation.
##If the user wants further knowledge about the simulation, please consider the following references.
##
##[1] T. Nakano and J.-Q. Liu, �Design and analysis of molecular relay
##channels: An information theoretic approach,� IEEE Transactions on
##NanoBioscience, vol. 9, pp. 213�221, 2010.
##
##[2] A. Goldbeter, G. Dupont, and M. J. Berridge, �Minimal model for signalinduced
##Ca2+ oscillations and for their frequency encoding through
##protein phosphorylation,� Proceedings of the National Academy of
##Science USA, vol. 87, pp. 1461�1465, 1990.
##
##[3] D. T. Gillespie, �Exact stochastic simulation of coupled chemical
##reactions,� Journal of Physical Chemistry, vol. 81, no. 25, pp. 2340�
##2361, 1977.
##
##[4] M. T. Barros, S. Balasubramaniam, B. Jennings, and Y. Koucheryavy,
##�Transmission protocols for calcium signaling based molecular communications
##in deformable cellular tissues,� Manuscript sumbitted for
##journal publication, 2014.
##
##[5] M. T. Barros, S. Balasubramaniam, and B. Jennings, �Error control for
##calcium signaling based molecular communications,� in 47th Annual
##Asilomar Conference on Signals, Systems, and Computers, 2013.
##
##[6] M. T. Barros, S. Balasubramaniam, and B. Jennings, �Integrating information
##metrics and molecular communications for inferencing and
##detecting cellular tissue deformation,� Manuscript sumbitted for journal
##publication, 2014.
##
###########################################################################################################################################
###########################################################################################################################################

from numpy import fabs
#import stochpy
#from cell import cell
from random import uniform
##from pylab import figure,plot,show
from numpy.random import multinomial, exponential,random
from numpy import arange, array, empty,zeros,log
##from pykov import *
##import gc
##gc.disable()

##################################################################################
##Class responsible for the guillespie method implementation
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

				#pv = abs(array([eq() for eq in pvi]))# #propensity vector
				a0 = pv.sum() #sum of all transition probabilities
		#                print tim, pv, a0
				tau = (-1/a0)*log(random())
				#tau = (-1/a0)*log(uniform(0.0,1.0))
				p='n'
				while p == 'n':
					##event = rout2(pv,a0) inserted the roulet function here 21/03/2014
					###############################
					#####roulete function
					##############################
					limite = a0*random()
					#print limite
					aux = 0
					k = 0
					while((k<len(pv)) and (aux<limite)):
								aux = aux + pv[k]
								k = k + 1
					k = k - 1
					event = k
					#event = multinomial(1,pv/a0) # event which will happen on this iteration
					#print event
					p2 = uniform(0.0,1.0)
					##p2 = random()
					sum1=0
					sum2=0
					for i in xrange(event):
						if l == i:
							sum1=sum2
						sum2+=pv[i]/a0
					if sum1 < p2 and sum2 >= p2:
						p='s'
				#ini += tm[:,event.nonzero()[0][0]]
				#print tc, ini
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
        #self.steps=steps

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

        #pv = abs(array([eq() for eq in pvi]))# #propensity vector
        a0 = pv.sum() #sum of all transition probabilities
#                print tim, pv, a0
        tau = (-1/a0)*log(random())
        #tau = (-1/a0)*log(uniform(0.0,1.0))
        p='n'
        while p == 'n':
            event = self.roulette(pv,a0)
            #event = multinomial(1,pv/a0) # event which will happen on this iteration
            #print event
            p2 = uniform(0.0,1.0)
            sum1=0
            sum2=0
            for i in xrange(event):
                if l == i:
                    sum1=sum2
                sum2+=pv[i]/a0
            if sum1 < p2 and sum2 >= p2:
                p='s'
        #ini += tm[:,event.nonzero()[0][0]]
        #print tc, ini
        tc += tau
        if self.debug == 's':
            print 'ini :',ini
            print 'pv: ',pv
            print 'R :',event
            print 'Value: ', pv[event]
            print 'Release: ', pv[400]
            print 'RXACT: ', pv[416]
            print 'RXINACT', pv[417]
        #print tau
        #steps +=1
        #self.res[tim,:,round] = ini
#        tvec = tvec[:tim]
#        self.res = res[:tim,:,round]
        return tau, event, pv

    def roulette(self, popsel, popnumber):
        limite = popnumber*random()
        #print limite
        aux = 0
        k = 0
        while((k<len(popsel)) and (aux<limite)):
                    aux = aux + popsel[k]
                    k = k + 1
        k = k - 1
        #print k
        #print len(popsel)
        return k

##################################################################################
#Main body of the simulator
##################################################################################

def main(destination,conc,n,m,x,txconc,deltaamp):
    print '###################################################'
    print 'Simulation of: 1) Destination: ', destination
    print '2) Conc: ', txconc, '3) 3D Tissue Size: ',n, m
    print '4) Length of the cell :',x
    print '###################################################'

    ##################################################################################
    ##Variable Initialization
    ##################################################################################

    #### astrocytes

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

    ####ephitelium

    #phh = [0.333333333333,0.397282861581,0.455089728717,0.507341806804,0.554575126347,0.597271531659,0.635867778391,0.670757405055,0.702295869972
    #,0.730805050486,0.756576081977,0.779871798426,0.800929736755,0.819964821146,0.837171500117,0.852725517033,0.866785620833,0.879495252417
    #,0.890984094151,0.901369413764,0.910757226704,0.919243350342,0.926914384574,0.933848611282,0.940116800537,0.945782928936,0.950904826691
    #,0.955534766602,0.95971999984,0.963503240068,0.966923099248,0.970014480683,0.972808934745,0.975334981116,0.977618400248,0.979682496561
    #,0.981548336139,0.983234961603,0.984759586465,0.986137770873,0.98738358045,0.988509729834,0.989527712414,0.990447917624,0.991279736989
    #,0.99203165998,0.992711360657,0.993325775987,0.993881176636,0.994383230972]

    #plh = [0.333333333333,0.301358569209,0.272455135642,0.246329096598,0.222712436827,0.201364234171,0.182066110804,0.164621297473,0.148852065014
    #,0.134597474757,0.121711959011,0.110064100787,0.0995351316227,0.090017589427,0.0814142499414,0.0736372414836,0.0666071895837,0.0602523737914
    #,0.0545079529247,0.0493152931181,0.044621386648,0.0403783248289,0.0365428077129,0.0330756943588,0.0299415997315,0.0271085355318,0.0245475866544
    #,0.022232616699,0.02014000008,0.0182483799658,0.016538450376,0.0149927596585,0.0135955326276,0.012332509442,0.0111907998759,0.0101587517193
    #,0.00922583193032,0.00838251919825,0.00762020676744,0.00693111456375,0.00630820977512,0.00574513508305,0.00523614379316,0.00477604118799
    #,0.0043601315054,0.00398417000987,0.00364431967125,0.00333711200656,0.00305941168212,0.00280838451412]

    #phl = [0.333333333333,0.301358569209,0.272455135642,0.246329096598,0.222712436827,0.201364234171,0.182066110804,0.164621297473,0.148852065014
    #,0.134597474757,0.121711959011,0.110064100787,0.0995351316227,0.090017589427,0.0814142499414,0.0736372414836,0.0666071895837,0.0602523737914
    #,0.0545079529247,0.0493152931181,0.044621386648,0.0403783248289,0.0365428077129,0.0330756943588,0.0299415997315,0.0271085355318,0.0245475866544
    #,0.022232616699,0.02014000008,0.0182483799658,0.016538450376,0.0149927596585,0.0135955326276,0.012332509442,0.0111907998759,0.0101587517193
    #,0.00922583193032,0.00838251919825,0.00762020676744,0.00693111456375,0.00630820977512,0.00574513508305,0.00523614379316,0.00477604118799
    #,0.0043601315054,0.00398417000987,0.00364431967125,0.00333711200656,0.00305941168212,0.00280838451412]

    #### xenopus

    #phh = [0.333333333333,0.410670791897,0.445466481454,0.459949002708,0.465060955179,0.466095442448,0.465539407907,0.464510810089,0.463476262365
    #,0.462598451758,0.461912367819,0.46139993636,0.461026806327,0.460760941063,0.460574480329,0.460443258373,0.460351936427,0.460290173831,0.460248091825
    #,0.460217834804,0.460196242295,0.460182442727,0.460173543649,0.460166466152,0.460161435794,0.460159164029,0.460157778359,0.460155703995
    #,0.460154290225,0.460154493019,0.460154435031,0.460153236783,0.460152782593,0.460153631556,0.460153681023,0.460152600557,0.460152545235
    #,0.460153562964,0.460153431072,0.46015237422,0.46015259763,0.460153606747,0.460153258821,0.460152304016,0.460152744974,0.460153597697
    #,0.460153073675,0.460152355963,0.460152931904,0.460153565747]

    #plh=[0.333333333333,0.361089611014,0.397501501773,0.43132158255,0.459003096423,0.480247650751,0.495947306903,0.507273290187,0.515312021
    #,0.520954695319,0.524882908321,0.527601928935,0.529476762612,0.530764843121,0.531647371462,0.532252137494,0.532665266546,0.53294618489
    #,0.533138132293,0.533270743029,0.533361828898,0.533422789776,0.533463623618,0.533492425246,0.53351229424,0.533524775959,0.533533137252
    #,0.533539997497,0.533544675711,0.533546751083,0.533548360069,0.533550616284,0.533551800574,0.533551482084,0.533551784093,0.533553089103
    #,0.533553310328,0.533552428664,0.533552637551,0.533553729952,0.533553548752,0.533552585366,0.533552944722,0.533553893926,0.533553470014
    #,0.533552639554,0.533553157692,0.533553864552,0.533553301572,0.533552681178]

    #phl=[0.333333333333,0.22823959709,0.157032016773,0.108729414741,0.0759359483979,0.0536569068007,0.0385132851899,0.0282158997238,0.021211716635
    #,0.0164468529234,0.0132047238598,0.010998134705,0.00949643106087,0.00847421581611,0.00777814820948,0.00730460413251,0.00698279702694
    #,0.00676364127942,0.00661377588202,0.00651142216702,0.0064419288067,0.00639476749629,0.00636283273226,0.00634110860156,0.00632626996593
    #,0.00631606001139,0.00630908438871,0.00630429850771,0.00630103406485,0.0062987558978,0.00629720490004,0.00629614693332,0.00629541683331
    #,0.00629488635978,0.00629453488348,0.00629431034042,0.00629414443703,0.0062940083719,0.00629393137696,0.00629389582846,0.00629385361765
    #,0.00629380788715,0.00629379645689,0.00629380205829,0.00629378501136,0.00629376274815,0.00629376863346,0.00629377948584,0.00629376652461
    #,0.0062937530747]

    D = 122500
    l = 0.5
    space = 10
    #laticeset = [cell()] * int(space/l)
    laticeset = [0] * int(space/l) * n * m
    t=0
    tEnd=5
    tPulse = 5
    tini=0
    tend=5
    C = []
    C2 = []
    C4 = []
    T = []
    #W=0.5#
    numberofrec = 29
    numberofvar = 30
    #ini = [100,10,1000,50000,50,2,1000,500000,50,2,4,2000,90,1,4.2,250,4,500,0.85,50,100000,0.5,0.6] * int(space/l)
    ini = [100,0,0.05,15,0.1,2.02,0.1,4000,50,1.5,50,0.3,8,0.05,0.15,0.15,0.1,0.1,4,0,0,0.5,D,(3.1416*((float(x)/2.0)**2)),0.0006,2.5,2.2,phh[0],phl[0],plh[0]] * int(space/l) * n * m
    vsars = ['0v1','1Y','2vin','3VM2','4C','5n','6K2','7VM3','8kout','9S','10kf','11kp','12kdeg','13vp','14kcaaa','15kcai','16kip3','17Z'
    ,'18q','19W','20A','21kia','22D','23l','24K','25ka','26m','27phh','28phl','29plh'] * int(space/l) * n * m
    prop = []
    freq = [0]*(numberofrec)
    freq2 = [0]*int(space/l) * n * m
    #ALPHA = 1/(6.02*100*(3.1416*(float(x)/2.0)**2)**3)
    ALPHA = 0.01
    q = 4
    freqi = []
    f=0
    f2=0
    re=0.2
    #logfile = open ( 'log.txt', 'w' )
    states=0
    stater=1
    liststates=[]
    liststater=[]
    ini[(len(laticeset)/2)*numberofvar+1]=txconc
    ini[(len(laticeset)/2)*numberofvar+4]=deltaamp
    sourceid = (len(laticeset)/2)*numberofvar+4
    dest = destination
    ini[(len(laticeset)/2+dest)*numberofvar+19]=conc
    destinationid = (len(laticeset)/2+dest)*numberofvar+4
    listdiffTx = [0]*4
    listdiffRx = [0]*4
    #log='n'
    debug='n'
    p=[]
    p2=[]
    p4=[]
    p5=[]
    lpn=[]
    lpn2=[]
    bit=[]
    error=[]
    Tb = 5
    concs = 0
    bits = []
    concscount=0
    f2=0
    #detecthers=(5.5*txconc)+0.045
    detecthers=0.2
    maxconc = 10000
    concdelay = 0
    maxGainSamples = 10
    GainTrans = []
    GainRec = []

    delay = 0
    totalcalciumtrans=0
    totalcalcium=deltaamp

    #### noise variables ###

    sourcenoise = 0
    destinationnoise = 0
    systemnoise = 0
    reflectionnoise = 0

    comsourcenoise = []
    comdestinationnoise = []
    comsystemnoise = []
    comreflectionnoise = []

    #logc2 = open ( 'logc'+ str(destination)+str(txconc) +'.txt', 'w' )

    ### light of sight list ###

    loslist = []
    loslistappend = loslist.append

    for c in xrange((len(laticeset)/2)+1,(len(laticeset)/2)+dest):
        loslistappend(c)

    print loslist


    a = 1.05 # linear coeficient

    props1=0
    props0=0
    propr1=0
    propr0=0

    probr1s0=0
    probr0s0=0
    probr1s1=0
    probr0s1=0

    probtest=0.5
    propappend = prop.append
    propextend = prop.extend

    ##################################################################################
    ## computing the length with pressing
    ##################################################################################

##    pressrate=0.3
##    pressarea=80
##    Area=x*x
##    a=(x-pressrate)
##
##    increrate=(x-(x-pressrate))/pressarea
##
##    for lol in range((len(ini)/2),pressarea*23+((len(ini)/n)*(n/2)),23):
##        ini[lol+21]=(2*Area)/(a+a+increrate)
##        ini[(len(ini)/2-23)-(lol-(len(ini)/2))+21]=(2*Area)/(a+a+increrate)
##        ini[lol+(23 * int(space/l))+21] = (2*Area)/(a+a+increrate)
##        ini[(len(ini)/2-23)-(lol-(len(ini)/2))+(23 * int(space/l))+21]=(2*Area)/(a+a+increrate)
##        ini[lol-(23 * int(space/l))+21] = (2*Area)/(a+a+increrate)
##        ini[(len(ini)/2-23)-(lol-(len(ini)/2))-(23 * int(space/l))+21]=(2*Area)/(a+a+increrate)
##        a=a+increrate

####    for cu in range(0,40):
####        if (cu>20 and cu<60) or (cu>100 and cu<140) or (cu>180 and cu<220):
####            ini[cu*23+21]=0.8


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
            print '1'
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
            propappend(lambda ini,i:(((ini[(i/numberofrec)*numberofvar+22]/(ini[(i/numberofrec)*numberofvar+23]))*fabs(ini[(i/numberofrec)*numberofvar+4]-ini[(i/numberofrec+(int(space/l)*n))*numberofvar+4]))*ini[(i/numberofrec)*numberofvar+28]) if ini[(i/numberofrec)*numberofvar+4] > ini[(i/numberofrec-int(space/l))*numberofvar+4] else 0)
        else:
            propappend(lambda ini,i:0)

        #plh
        if tant+(numberofvar * int(space/l) * n) < len(ini):
            propappend(lambda ini,i:(((ini[(i/numberofrec)*numberofvar+22]/(ini[(i/numberofrec)*numberofvar+23]))*fabs(ini[(i/numberofrec)*numberofvar+4]-ini[(i/numberofrec+(int(space/l)*n))*numberofvar+4]))*ini[(i/numberofrec)*numberofvar+29]) if ini[(i/numberofrec)*numberofvar+4] > ini[(i/numberofrec-int(space/l))*numberofvar+4] else 0)
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

    ##logfilewrite = logfile.write
    while t <= tEnd:

        r = [.001,.1]
        tm = 2

        M = Model(vnames = vars,rates = r,inits=ini, tmat=tm,propensity=prop,debug=debug)
        tau,event,pv = M.run(tmax=200,reps=1)
        #print t
        t = t + tau*1000
        #reaction = event.nonzero()[0][0]
        reaction = event
        #print reaction
        i=reaction/numberofrec
        #print i
        #print reaction
        x=(reaction-(reaction/numberofrec)*numberofrec)+1

##        if reaction>100:
##            x = int(str(reaction+1)[2])
##        if reaction<100 and reaction>9:
##            x = int(str(reaction+1)[1])
##        if reaction<10:
##            x = int(str(reaction+1))

        #print x
        #print vsars[i*23+4]

        #if log == 's':
            #logfilewrite('Time : '+str(t))
            #logfilewrite(' Reaction : '+str(x))
            #logfilewrite('Cell : '+str(i))

        if int(t) != f:
            f = int(t)
            if int(t) >= tini and int(t) <= tend:
                ini[(len(laticeset)/2)*numberofvar+4]=deltaamp
            print 'Time:',f
            print 'A :', ini[(len(laticeset)/2+dest)*numberofvar+20]
            print 'Freq: ',freq
            #print 'Reaction in the cell: ',freq2
            print '###### numero difusoes'
            print listdiffTx
            print listdiffRx
            print '######### Noises ####'
            print 'Source Noise: ',sourcenoise
            print 'Destination Noise: ',destinationnoise
            print 'System Noise: ',systemnoise
            print 'Reflection Noise: ',reflectionnoise
            print '#####################'
            print totalcalcium
            print sum(C)
            print tEnd
            
            print '#######new method#####'
            print deltaamp
            print concdelay

            #if int(t) == 2 or int(t) == 4 or int(t) == 6 or int(t) == 10:
            if int(t) > 0:
				a= 0
				for j in range(0,m):
					for k in range(0,int(space/l)*numberofvar,numberofvar):
						for i in range(0,n):
							print ini[k+((a+i)*int(space/l)*numberofvar)+4],
						print ''
					a=a+n
					print ''

            itme = 0
            comsourcenoise.append(sourcenoise)
            comdestinationnoise.append(destinationnoise)
            comsystemnoise.append(systemnoise)
            comreflectionnoise.append(reflectionnoise)
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

			##bit detection
            if f%Tb == 0  and f!=0:
                mean=float(concs)/float(concscount)
                if mean > detecthers:
                    bits.append('1')
                else:
                    bits.append('0')
                concs=0
                concscount=0
            #if t >10 and t < 20 :
            #    ini[(len(laticeset)/2)*23+4]+=ini[(len(laticeset)/2)*23+1]*ini[(len(laticeset)/2)*23]

        #if x == 1 and ((t >= 10 and t <= 20) or (t >= 30 and t <= 40) or (t >= 50 and t <= 60) or (t >= 70 and t <= 80) or (t >= 90 and t <= 100) or (t >= 110 and t <= 120) or (t >= 130 and t <= 140) or (t >= 150 and t <= 160) or (t >= 170 and t <= 180) or (t >= 190 and t <= 200)):
        #if x == 1 and t >10 and t < 20 and re < ini[0]*ini[1]:
        #if x == 1 and (t >= tini and t <= tend) and totalcalcium <= maxconc:
            ##ini[i*23+4]+=ALPHA
            #re+=ALPHA
            ##if int(t) != f2:
            ##    f2 = int(t)
            #ini[(len(laticeset)/2)*numberofvar+4]+=ini[(len(laticeset)/2)*numberofvar+1]*ini[(len(laticeset)/2)*numberofvar]-ini[(len(laticeset)/2)*numberofvar+4]
            #totalcalcium=totalcalcium+ini[(len(laticeset)/2)*numberofvar+4]
            ##ini[i*23+4]+=ALPHA


        if t >= tini and t <= tend+1:
            if ini[(len(laticeset)/2+dest)*numberofvar+4] > detecthers:
                bit.append('1')
            else:
                bit.append('0')
                error.append('1')
        else:
            if ini[(len(laticeset)/2+dest)*numberofvar+4] < detecthers:
                bit.append('0')
            else:
                 bit.append('1')
                 error.append('1')

        concs = concs +  ini[(len(laticeset)/2+dest)*numberofvar+4]
        concscount = concscount+1

        #if int(t) == 10 or int(t) == 30 or int(t) == 50 or int(t) == 70 or int(t) == 90 or int(t) == 110 or int(t) == 130 or int(t) == 150 or int(t) == 170 or int(t) == 190:
        if int(t) == tini:
            states=1
        #if int (t) == 21 or int (t) == 41 or int (t) == 61 or int (t) == 81 or int (t) == 101 or int (t) == 121 or int (t) == 141 or int (t) == 161 or int (t) == 181:
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

        #if x == 7 and stater == 0:
        if x == 10 :
            ini[i*numberofvar+4]-=q*ALPHA
            ini[i*numberofvar+19]+=q*ALPHA
            if i == len(laticeset)/2+dest:
                stater=1
                ini[(len(laticeset)/2+dest)*numberofvar+20]+=1
            if ini[i*numberofvar+4] < 0:
                ini[i*numberofvar+4] = ini[i*numberofvar+4]*(-1)
            #print t
        #if x == 8 and stater == 11
        if x == 11 :
            ini[i*numberofvar+4]+=q*ALPHA
            ini[i*numberofvar+19]-=q*ALPHA
            if i == len(laticeset)/2+dest:
                stater=0
                ini[(len(laticeset)/2+dest)*numberofvar+20]-=1
            if ini[i*numberofvar+19] < 0:
                ini[i*numberofvar+19] = ini[i*numberofvar+19]*(-1)
            #print t91

        ## DIFFUSION REACTIONS

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

            ## calculo do delay
            if ((i+1)*numberofvar+4) == destinationid:
				concdelay +=ALPHA 
            
            ## source noise ##

            if (i < len(laticeset)-1) and (((i+1)*numberofvar+4) == sourceid) :
                sourcenoise = sourcenoise + ALPHA

            ## destination noise ##soilwo

            if (i != 0) and (((i-1)*numberofvar+4) == destinationid) and (t < tini or t > tend) :
                destinationnoise = destinationnoise + ALPHA

            if ((i*numberofvar+4) == destinationid) :
                destinationnoise = destinationnoise + ALPHA

            ## system noise  ##

            if ((i-1) in loslist) and (t < tini or t > tend) :
                systemnoise = systemnoise+ALPHA

            if ((i) in loslist) and (t < tend or t > tini) :
                systemnoise = systemnoise+ALPHA

            ## reflection noise ##


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

            ## calculo do delay
            if ((i-1)*numberofvar+4) == destinationid:
				concdelay +=ALPHA 
            
            ## source noise ##

            if (i*numberofvar+4) == sourceid :
                sourcenoise = sourcenoise + ALPHA

            ## destination noise ##

            if (i != 0) and (((i-1)*numberofvar+4) == destinationid) :
                destinationnoise = destinationnoise + ALPHA

            if ((i*numberofvar+4) == destinationid) and (t < tini or t > tend) :
                destinationnoise = destinationnoise + ALPHA

            ## system noise  ##

            if ((i-1) in loslist) and (t < tini or t >tend) :
                systemnoise = systemnoise+ALPHA

            if ((i) in loslist) and (t < tini or t > tend) :
                systemnoise = systemnoise+ALPHA

            ## reflection noise ##

            if ((i) in loslist) and (t >= tini and t <= tend) :
                reflectionnoise = reflectionnoise+ALPHA


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

            ## calculo do delay
            if ((i+int(space/l))*numberofvar+4) == destinationid:
				concdelay +=ALPHA 
            
            ## source noise ##

            if (i*numberofvar+4) == sourceid :
                sourcenoise = sourcenoise + ALPHA

            if (i < len(laticeset)-1) and (((i+int(space/l))*numberofvar+4) == sourceid) :
                sourcenoise = sourcenoise + ALPHA

            ## destination noise ##

            if (i*numberofvar+4) == destinationid :
                destinationnoise = destinationnoise + ALPHA

            if (i < len(laticeset)-1) and (((i+int(space/l))*numberofvar+4) == destinationid) :
                destinationnoise = destinationnoise + ALPHA

            ## system noise  ##

            if ((i+int(space/l)) in loslist) and (t < tini or t > tend) :
                systemnoise = systemnoise+ALPHA

            if ((i) in loslist) and (t < tini or t > tend) :
                systemnoise = systemnoise+ALPHA

            ## reflection noise ##

            if ((i) in loslist) and (t >= tini and t <= tend) :
                reflectionnoise = reflectionnoise+ALPHA


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

            ## calculo do delay
            if ((i-int(space/l))*numberofvar+4) == destinationid:
				concdelay +=ALPHA 
            
            ## source noise ##

            if (i*numberofvar+4) == sourceid :
                sourcenoise = sourcenoise + ALPHA

            if (i != 0) and (((i-int(space/l))*numberofvar+4) == sourceid) :
                sourcenoise = sourcenoise + ALPHA

            ## destination noise ##

            if (i*numberofvar+4) == destinationid :
                destinationnoise = destinationnoise + ALPHA

            if (i != 0) and (((i-int(space/l))*numberofvar+4) == destinationid) :
                destinationnoise = destinationnoise + ALPHA

            ## system noise  ##

            if ((i-(int(space/l)*n)) in loslist) and (t < tini or t > tend) :
                systemnoise = systemnoise+ALPHA

            if ((i) in loslist) and (t < tini or t > tend) :
                systemnoise = systemnoise+ALPHA

            ## reflection noise ##

            if ((i) in loslist) and (t >= 10 and t <= 20) :
                reflectionnoise = reflectionnoise+ALPHA

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

            ## calculo do delay
            if ((i+(int(space/l)*n))*numberofvar+4) == destinationid:
				concdelay +=ALPHA 
            
            ## source noise ##

            if (i*numberofvar+4) == sourceid :
                sourcenoise = sourcenoise + ALPHA

            if (i < (int(space/l)*n)-1) and (((i+(int(space/l)*n))*numberofvar+4) == sourceid) :
                sourcenoise = sourcenoise + ALPHA

            ## destination noise ##

            if (i*numberofvar+4) == destinationid :
                destinationnoise = destinationnoise + ALPHA

            if (i < len(laticeset)-1) and (((i+(int(space/l)*n))*numberofvar+4) == destinationid) :
                destinationnoise = destinationnoise + ALPHA

            ## system noise  ##

            if ((i+int(space/l)) in loslist) and (t < tini or t > tend) :
                systemnoise = systemnoise+ALPHA

            if ((i) in loslist) and (t < tini or t > tend) :
                systemnoise = systemnoise+ALPHA

            ## reflection noise ##

            if ((i) in loslist) and (t >= tini and t <= tend) :
                reflectionnoise = reflectionnoise+ALPHA


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

            ## calculo do delay
            if ((i-(int(space/l)*n))*numberofvar+4) == destinationid:
				concdelay +=ALPHA 
            
            ## source noise ##

            if (i*numberofvar+4) == sourceid :
                sourcenoise = sourcenoise + ALPHA

            if (i != 0) and (((i-(int(space/l)*n))*numberofvar+4) == sourceid) :
                sourcenoise = sourcenoise + ALPHA

            ## destination noise ##

            if (i*numberofvar+4) == destinationid :
                destinationnoise = destinationnoise + ALPHA

            if (i != 0) and (((i-(int(space/l)*n))*numberofvar+4) == destinationid) :
                destinationnoise = destinationnoise + ALPHA

            ## system noise  ##

            if ((i-int(space/l)) in loslist) and (t < tini or t > tend) :
                systemnoise = systemnoise+ALPHA

            if ((i) in loslist) and (t < tini or t > tend) :
                systemnoise = systemnoise+ALPHA

            ## reflection noise ##

            if ((i) in loslist) and (t >= 10 and t <= 20) :
                reflectionnoise = reflectionnoise+ALPHA


        freq[x-1]= freq[x-1] + 1
        #freqi.append(i)

        #if log == 's':
            #logfile.write(' C state : '+str(ini[(len(laticeset)/2+dest)*numberofvar+4]))
            ##logfile.write( Ini : '+str(ini))
            #logfile.write('\n')

        #logc2.write(str(ini[(len(laticeset)/2+destination)*numberofvar+4])+'\n')

        C.append(ini[(len(laticeset)/2+dest)*numberofvar+4])
        C2.append(ini[(len(laticeset)/2)*numberofvar+4])
##        C4.append(ini[(len(laticeset)/2+9)*23+4])
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

        if (float(liststater.count(1))/float(len(liststater))*(ini[(len(laticeset)/2+dest)*numberofvar+4])) != 0:
            p4.append( ( (ini[(len(laticeset)/2+dest)*numberofvar+4]) / (ini[(len(laticeset)/2)*numberofvar+1]*ini[(len(laticeset)/2)*numberofvar])) )
        else:
            p4.append(0)

        #props1=float(sumpx1)/float(len(px1))
        props1=float(sum(liststates))/float(len(liststates))
        props0=1-props1
        propr1=float(liststater.count(1))/float(len(liststater))
        propr0=1-propr1

        if propr0 != 0 and propr1 != 0:
            probr1s0=(props0*propr1)/propr1
            probr0s0=1-probr1s0
            probr1s1=(props1*propr1)/propr1
            #probr1s1=propr1+probr1s0-(propr1*probr1s0)
            probr0s1=1-probr1s1


        p5.append(propr1*probr1s1)

        if p4[-1] >= 0.2 :
                     lpn.append(1)
        else:
			lpn.append(0)
        ##################################################################################
        ## delay
        ##################################################################################

        #if delay == 0 and ini[(len(laticeset)/2+dest)*numberofvar+4] <= 0.3 and tend <= t:
        #    delay = t

        #if t < tPulse :
        #    totalcalcium = sum(C2)

        #if concdelay > (float(deltaamp)/10.0):
        if sum(C) >= sum(C2):
            delay = t
            tEnd = 0
        else:
            tEnd+=tau*1000

    print '###### Number of Diffusions'
    print listdiffTx
    print listdiffRx
    print '#########'

    print '##### Frequencies'
    print freq
    print '#########'

##    print ""
##    print "###### METODO PROB 4 #####"
##    print ""
##
##    px1 = float(sum(liststates))/float(len(liststates))
##    px0 = 1 - px1
##
##    py1 = float(sum(p4))/float(len(p4))
##    py0 = 1 - py1
##
##    py1x0 = float(sum(p2))/float(len(p2))
##    py0x0 = 1 - py1x0
##
##    if py1x0==0:
##        py1x0=0.00000000001
##        py0x0 =0.9999999999
##
##    py1x1 = float(sum(p))/float(len(p))
##    py0x1 = 1 - py1x1
##
##    print px0,px1,py0,py1
##    print py1x0,py0x0,py1x1,py0x1
##    print py1*py1x0,py1*py0x0,py1*py1x1,py1*py0x1
##
##    from math import log
##
##    ibits= px0*py0x0*(log(py0x0/(py0))/log(2))+px1*py0x1*(log(py0x1/(py0))/log(2))+px0*py1x0*(log(py1x0/(py1))/log(2))+px1*py1x1*(log(py1x1/(py1))/log(2))
##
##    print ibits

    print ""
    print "###### METHOD PROB 42 #####"
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

    #py1x1 = (px1*py1)/py1
    py1x1 = py1+py1x0-(py1*py1x0)
    py0x1 = 1 - py1x1

    print px0,px1,py0,py1
    print py1x0,py0x0,py1x1,py0x1
    print py1*py1x0,py1*py0x0,py1*py1x1,py1*py0x1

    from math import log

    ibits= px0*py0x0*(log(py0x0/(py0))/log(2))+px1*py0x1*(log(py0x1/(py0))/log(2))+px0*py1x0*(log(py1x0/(py1))/log(2))+px1*py1x1*(log(py1x1/(py1))/log(2))

    print ibits

    print 'P[error] = ', (float(len(error))/float(len(bit)))
    print len(error)
    print len(bit)
    print bits
    print 'Concentration', float(sum(C))/float(len(C))
    print 'Gain', 10*log((float(sum(C))/float(len(C)))/float((float(sum(C2))/float(len(C2)))))
    print 'Gain', 10*log(max(C))/(max(C2))
    print 'Delay', t
    print 'newGAIN', 10*log((float(sum(GainRec))/float(len(GainRec)))/(float(sum(GainTrans))/float(len(GainTrans))) )

    print '######### Noises ####'
    print 'Source Noise: ',sourcenoise
    print 'Destination Noise: ',destinationnoise
    print 'System Noise: ',systemnoise
    print 'Reflection Noise: ',reflectionnoise
    print '#####################'

    print '######### comulative Noises ####'
    print 'Source Noise: ',comsourcenoise
    print 'Destination Noise: ',comdestinationnoise
    print 'System Noise: ',comsystemnoise
    print 'Reflection Noise: ',comreflectionnoise
    print '#####################'

    # transfile = open ( 'Trans-'+str(destination)+str(conc)+str(n)+str(m)+str(x)+str(txconc)+str(deltaamp)+'.txt', 'w' )
    # recfile = open ( 'Rec-'+str(destination)+str(conc)+str(n)+str(m)+str(x)+str(txconc)+str(deltaamp)+'.txt', 'w' )
    
    # transfile.write('[')
    # recfile.write('[')
    
    # print ''
    # print '################ Concs#####'
    # print 'Trans'
    # for i in C2:
    #     transfile.write(str(i)+',')
    # print ''
    # print 'Receiv'
    # for j in C:
    #     recfile.write(str(j)+',')
    
    # transfile.write(']')
    # recfile.write(']')
    
    # transfile.close()
    # recfile.close()

    #logc2.close()

    #logfile.close()
    #figure(1)
    #plot(T,C)
    #figure(2)
    #plot(T,C2)
####    figure(3)
####    plot(T,C4)
####    figure(4)
####    plot(T,p4)
####    figure(5)
####    plot(T,p5)
##    savefig('graphTxconc'+str(txconc)+'.pdf')
##    show()

    ##print bit

    ##print liststates
    ##print liststater

##    ##
##    ## Gilbert-Elliot Model
##    ##
##
##    ## We represent the states of x and y with the following: X0Y0 -> a, X1Y1 -> b, X0Y1 -> c, X1Y0 -> d. The probabilities of each state are pi, where i is the state letter.
##    ## The transition probabilities are pij, where i is the current state and y is the next state
##
##    pa=0
##    pb=0
##    pc=0
##    pd=0
##
##    paa=0
##    pab=0
##    pac=0
##    pad=0
##
##    pba=0
##    pbb=0
##    pbc=0
##    pbd=0
##
##    pca=0
##    pcb=0
##    pcc=0
##    pcd=0
##
##    pda=0
##    pdb=0
##    pdc=0
##    pdd=0
##
##    ## calculate the probabilities of each state
##
##    listpa = 0
##    listpb = 0
##    listpc = 0
##    listpd = 0
##
##    for i in range(0,len(liststates)):
##        if liststates[i] == 0 and liststater[i] == 0 :
##            listpa = listpa + 1
##        if liststates[i] == 1 and liststater[i] == 1 :
##            listpb = listpb + 1
##        if liststates[i] == 0 and liststater[i] == 1 :
##            listpc = listpc + 1
##        if liststates[i] == 1 and liststater[i] == 0 :
##            listpd = listpd + 1
##
##    print len(liststates)
##    print len(liststater)
##
##    pa = float(listpa)/float(len(liststates))
##    pb = float(listpb)/float(len(liststates))
##    pc = float(listpc)/float(len(liststates))
##    pd = float(listpd)/float(len(liststates))
##
##    print ''
##    print '#########State Probabilities########'
##    print ''
##
##    print pa
##    print pb
##    print pc
##    print pd
##
##    ## calculate the transistion probabilities
##
##    listpaa=0
##    listpab=0
##    listpac=0
##    listpad=0
##
##    listpba=0
##    listpbb=0
##    listpbc=0
##    listpbd=0
##
##    listpca=0
##    listpcb=0
##    listpcc=0
##    listpcd=0
##
##    listpda=0
##    listpdb=0
##    listpdc=0
##    listpdd=0
##
##    for i in range(1,len(liststates)):
##        if liststates[i] == 0 and liststater[i] == 0 and liststates[i-1] == 0 and liststater[i-1] == 0 :
##            listpaa = listpaa + 1
##        if liststates[i] == 0 and liststater[i] == 0 and liststates[i-1] == 1 and liststater[i-1] == 1 :
##            listpab = listpab + 1
##        if liststates[i] == 0 and liststater[i] == 0 and liststates[i-1] == 0 and liststater[i-1] == 1 :
##            listpac = listpac + 1
##        if liststates[i] == 0 and liststater[i] == 0 and liststates[i-1] == 1 and liststater[i-1] == 0 :
##            listpad = listpad + 1
##
##        if liststates[i] == 1 and liststater[i] == 1 and liststates[i-1] == 0 and liststater[i-1] == 0 :
##            listpba = listpba + 1
##        if liststates[i] == 1 and liststater[i] == 1 and liststates[i-1] == 1 and liststater[i-1] == 1 :
##            listpbb = listpbb + 1
##        if liststates[i] == 1 and liststater[i] == 1 and liststates[i-1] == 0 and liststater[i-1] == 1 :
##            listpbc = listpbc + 1
##        if liststates[i] == 1 and liststater[i] == 1 and liststates[i-1] == 1 and liststater[i-1] == 0 :
##            listpbd = listpbd + 1
##
##        if liststates[i] == 0 and liststater[i] == 1 and liststates[i-1] == 0 and liststater[i-1] == 0 :
##            listpca = listpca + 1
##        if liststates[i] == 0 and liststater[i] == 1 and liststates[i-1] == 1 and liststater[i-1] == 1 :
##            listpcb = listpcb + 1
##        if liststates[i] == 0 and liststater[i] == 1 and liststates[i-1] == 0 and liststater[i-1] == 1 :
##            listpcc = listpcc + 1
##        if liststates[i] == 0 and liststater[i] == 1 and liststates[i-1] == 1 and liststater[i-1] == 0 :
##            listpcd = listpcd + 1
##
##        if liststates[i] == 1 and liststater[i] == 0 and liststates[i-1] == 0 and liststater[i-1] == 0 :
##            listpda = listpda + 1
##        if liststates[i] == 1 and liststater[i] == 0 and liststates[i-1] == 1 and liststater[i-1] == 1 :
##            listpdb = listpdb + 1
##        if liststates[i] == 1 and liststater[i] == 0 and liststates[i-1] == 0 and liststater[i-1] == 1 :
##            listpdc = listpdc + 1
##        if liststates[i] == 1 and liststater[i] == 0 and liststates[i-1] == 1 and liststater[i-1] == 0 :
##            listpdd = listpdd + 1
##
##    print ''
##    print '#########Transition Probabilities########'
##    print ''
##
##    paa=float(listpaa)/float(listpa)
##    pab=float(listpab)/float(listpa)
##    pac=float(listpac)/float(listpa)
##    pad=float(listpad)/float(listpa)
##
##    print paa
##    print pab
##    print pac
##    print pad
##
##    pba=float(listpba)/float(listpb)tEnd
##    pbb=float(listpbb)/float(listpb)
##    pbc=float(listpbc)/float(listpb)
##    pbd=float(listpbd)/float(listpb)
##
##    print pba
##    print pbb
##    print pbc
##    print pbd
##
##    pca=float(listpca)/float(listpc)
##    pcb=float(listpcb)/float(listpc)
##    pcc=float(listpcc)/float(listpc)
##    pcd=float(listpcd)/float(listpc)
##
##    print pca
##    print pcb
##    print pcc
##    print pcd
##
##    if listpd != 0:
##        pda=float(listpda)/float(listpd)
##        pdb=float(listpdb)/float(listpd)
##        pdc=float(listpdc)/float(listpd)
##        pdd=float(listpdd)/float(listpd)
##
##        print pda
##        print pdb
##        print pdc
##        print pdd
##    else:
##        print pda
##        print pdb
##        print pdc
##        print pdd
##
##    T = Chain({('A','A'): paa, ('A','B'): pab, ('A','C'): pac,('A','D'): pad,
##               ('B','A'): pba, ('B','B'): pbb,('B','C'): pbc, ('B','D'): pbd,
##               ('C','A'): pca,('C','B'): pcb, ('C','C'): pcc, ('C','D'): pcd,
##               ('D','A'): pda, ('D','B'): pdb, ('D','C'): pdc,('D','D'): pdd})
##
##    print T.entropy()
##    print T.entropy(norm=True)


##################################################################################
##Calling up the simulations
##################################################################################

main(1,10,3,3,5,0.1,1)
main(2,10,3,3,5,0.1,1)
main(3,10,3,3,5,0.1,1)
main(4,10,3,3,5,0.1,1)
main(5,10,3,3,5,0.1,1)
main(6,10,3,3,5,0.1,1)
main(7,10,3,3,5,0.1,1)
main(8,10,3,3,5,0.1,1)
main(9,10,3,3,5,0.1,1)

#main(1,0.5,3,0.5,0.01)
#main(13,0.5,3,3,0.5,0.01)
#main(3,0.5,3,3,0.5,0.01)
#main(5,0.5,3,3,0.5,0.01)
#main(9,0.5,3,3,0.5,0.01)
#main(13,0.5,3,3,0.5,0.01)

#main(13,0.5,3,3,5,0.)
#main(13,1,3,3,5,0.05)
#main(13,10,3,3,5,0.05)
#main(13,100,3,3,5,0.1)
#main(13,1000,3,3,5,0.05)
