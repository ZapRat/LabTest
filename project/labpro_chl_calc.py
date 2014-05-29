# -*- coding: utf-8 -*-
"""
Created on Sun Apr 20 16:05:53 2014
Script to run spectometry-based pigment concentrations.

@author: Charles
"""

import pdb
import json
import numpy as np
#general utilities
class TableMaker:
    def __init__(self,out_fn,header_str,format_str):
        """Start to write a table to a file named in out_fn.
        The first row will be a CSV string in header_str,
        each row will be formatted with format_str"""
        self.form=format_str
        self.writer=open(out_fn,'wt')
        self.writer.write(header_str+'\n')
    def writeRow(self,outline):
        """Write a line to the table; it should be a list with the
        same number of elements as template entries in the format string."""
        self.writer.write(self.form%tuple(outline)+"\n")
    def close(self):
        """Close the output file"""
        self.writer.close()
        
#spectrophotometer section
WAVELENGTHS=[750,664,647,630,510,480]
WAVES_AA=[750,664]#wavelengths used after acid
def getChlAspec(E,small_v,filtered_vol,cel_len=10):
    """Get the chlorophyll A based on the abosorption 
    spectrum represented by E, which is a dictionary
    where the wavelengths are the keys and the absorptions 
    are the values."""
    #last part of the lab equations concerning various 
    #volumes
    volfac=(10.0/cel_len)*(small_v/10.0)*(1.0/filtered_vol)
    chl_a=(11.85*E[664]-1.54*E[647]-.08*E[630])*volfac
    return chl_a
def getChlB(E,small_v,filtered_vol,cel_len=10):
    """Get the chlorophyll B based on the abosorption 
    spectrum represented by E, which is a dictionary
    where the wavelengths are the keys and the absorptions 
    are the values."""
    #last part of the lab equations concerning various 
    #volumes
    volfac=(10.0/cel_len)*(small_v/10.0)*(1.0/filtered_vol)
    return (-5.43*E[664]+21.03*E[647]-2.66*E[630])*volfac
    
def getChlC(E,small_v,filtered_vol,cel_len=10):
    """Get the chlorophyll C based on the abosorption 
    spectrum represented by E, which is a dictionary
    where the wavelengths are the keys and the absorptions 
    are the values."""
    #last part of the lab equations concerning various 
    #volumes
    volfac=(10.0/cel_len)*(small_v/10.0)*(1.0/filtered_vol)
    return (-1.67*E[664]-7.60*E[647]+24.52*E[630])*volfac

def getCarotenoids(E,small_v,filtered_vol,cel_len=10):
    """Get the caretenoids based on the spectrum."""
    volfac=(10.0/cel_len)*(small_v/10.0)*(1.0/filtered_vol)
    return 7.6*(E[480]-1.49*E[510])*volfac
    
def getChlAacid(E664b,E664a,small_v,filtered_vol,cel_len=10):
    """Get the Chl a using a technique that compares absorbance at 664
    before and after acidification"""
    return (26.7*(E664b-E664a)*small_v)/(filtered_vol*cel_len)

def getPhaeopig(E664b,E664a,small_v,filtered_vol,cel_len=10):
    """Get the Phaeopigments a using a technique that compares absorbance 
    at 664 before and after acidification"""
    return (26.7*(1.7*E664a-E664b)*small_v)/(filtered_vol*cel_len)
      
def getCorrectedAbs(E_in,blanks,cel_already=False):
    """Correct the spectrum E_in based on the blank spectrum, blanks. 
    If cel_already is 
    True, the correction only applies the turbidity blank, if it is False,
    Each value is corrected both by the blank for its frequency and turbidity"""
    #make the first round of correction based on blanks, if cel_already is false
    cel_corrected={}
    for curwave in WAVELENGTHS:
        if curwave in E_in:#after acid spectra are shorter, so skip unused 
            #wavelengths
            if cel_already:
                cel_corrected[curwave]=E_in[curwave]
            else:
                cel_corrected[curwave]=E_in[curwave]-blanks[curwave]
    #now do the turbidity correction
    #handle the 750 spectrum separately, it's special
    turbid_blank=cel_corrected[750]
    #initialize the fully corrected spectrum
    corrected_spec={750:cel_corrected[750]}
    #the top 3 below 750 get the turbidity blank subtracted from them
    for curwave in [664,647,630]:
        if curwave in E_in:#after acid spectra are shorter, so skip unused 
            corrected_spec[curwave]=cel_corrected[curwave]-turbid_blank
    #the shortes wavelenthse get 2x and 3x turbid blank correction
    if 510 in E_in:#after acid spectra are shorter, so skip unused 
        corrected_spec[510]=cel_corrected[510]-2*turbid_blank
    if 480 in E_in:
        corrected_spec[480]=cel_corrected[480]-3*turbid_blank
    return corrected_spec

def listToSpectrum(speclist):
    """Converts a list into a spectrum dictionary assuming that 
    the first item has a wavelength of 750, and each subsequent item 
    has a wavelength that matches it's index in the 'WAVELENGTHS' constant.
    speclist must be the same size or shorter than WAVELENGTHS, if it's 
    shorter, the lower wavelengths will be omitted from the output.
    The return value is a dictionary where the key is the wavelength
    and the value is the absorbance."""
    bcount=len(speclist)
    if bcount>len(WAVELENGTHS):
        raise Exception("More bands than standard wavelengths!")
    retriever={}#spec
    for i_band in range(bcount):
        retriever[WAVELENGTHS[i_band]]=speclist[i_band]
    return retriever

class PigSpectSite:
    """This class represents the spectrum data for a sampling 
    site having 3 replicates."""
    def __init__(self,sitenum,blanks):
        self.reps={}
        self.sitenum=sitenum
        self.blanks=blanks
        self.acetone_vol=50.0#all filteres were diluted to an aceton vol of
            #50mL

    def addRep(self,rep_index,rep_dict):
        self.reps[rep_index]=rep_dict
    def describe(self):
        """Dumps dat for this replicate to the screen"""
        rep_order=sorted(self.reps.keys())
        print "Site: %d================-----------"%self.sitenum
        for currep in rep_order:
            currep_dump=json.dumps(self.reps[currep],sort_keys=True,indent=4)
            print "Replicate %d:\n %s"%(currep,currep_dump)
            #pdb.set_trace()
            ba_dump=json.dumps(self.pigs_ba[currep],sort_keys=True,indent=4)
            print "Pigments from Before Acid Spectrum:\n%s"%ba_dump
            aa_dump=json.dumps(self.pigs_aa[currep],sort_keys=True,indent=4)
            print "Pigments from Before and After Acid Spectrum:\n%s"%aa_dump
            
        print "Averages--------------------------------"
        ba_avg_dump=json.dumps(self.avgpig_ba,
                                sort_keys=True,indent=4)
        print "Average Pigments from Before Acid Spectrum:\n%s"%ba_avg_dump
        aa_avg_dump=json.dumps(self.avgpig_aa,
                                sort_keys=True,indent=4)
        print "Average Pigments from Bef/Aft Acid Spectrum:\n%s"%aa_avg_dump

    def compute(self):
        """This computes all of the pigment concentrations and places them
        in local variables.  It should be called after all replicates have
        been added"""
        self.pigs_ba={}#this dictionary is keyed by replicate number
                #each element is a dictionary with each pigment name as 
                #the key and the concentration as the value.
                #these are the pigments computed by the spectrum
        self.pigs_aa={}#this dictionary is keyed by replicate number
                #each element is a dictionary with each pigment name as 
                #the key and the concentration as the value.
                #these are the pigments computed difference between the
                #spectra before and after acidification
        for repnum in sorted(self.reps.keys()):
            currdata=self.reps[repnum]
            
            corr_E=getCorrectedAbs(currdata['ba'],self.blanks)
            #the 'vol' member of the input is in mL, but the formulae need Liters,
            #so do a conversion
            vol_L=currdata['vol']*1e-3
            #do the before acid calculations and store them
            self.pigs_ba[repnum]={
            'chl_a':getChlAspec(corr_E,self.acetone_vol,vol_L),
            'chl_b':getChlB(corr_E,self.acetone_vol,vol_L),
            'chl_c':getChlC(corr_E,self.acetone_vol,vol_L)}

            corr_E_aa=getCorrectedAbs(currdata['aa'],self.blanks)
            self.pigs_aa[repnum]={
                'chl_a':getChlAacid(corr_E[664],corr_E_aa[664],
                                    self.acetone_vol,vol_L),    
                'phaeopig':getPhaeopig(corr_E[664],corr_E_aa[664],
                                       self.acetone_vol,vol_L)    
                
            }

        #now make averages of the replicates            
        self.avgpig_ba=avgReps(self.pigs_ba,['chl_a','chl_b','chl_c']) 
        self.avgpig_aa=avgReps(self.pigs_aa,['chl_a','phaeopig'])            

def avgReps(reps,pig_keys):
    """Average Replicate pigment data stored in dictionaries into
    a dictionary of the same structure.  pig_keys tells which keys 
    from the original dictionaries are to be averaged."""
    #now make averages of the replicates
    replists={}#will be lists of replicate value keyed by
        #which pigment concentration is being listed
    #make the lists
    for repnum in sorted(reps.keys()):
        for curpig in pig_keys:
            if curpig not in replists:#start the list to average
                replists[curpig]=[]
            replists[curpig].append(reps[repnum][curpig])
    #make the averages
    retriever={}
    for curpig in replists.keys():
        retriever[curpig]=np.mean(replists[curpig])
    return retriever
    
class PigBySpectrum:
    """This class processes all spectrophotometer-based pigment calculations"""        
    def getSpectrumInputs(self):
        """This routine contains the input data from before acidification; 
        it outputs the data in a way that's easy to use by the various pigment
        calculation routines""" 
        #the first element of the table is a tuple with the site and replicate
        spec_table=[
        [(1,1),-.03 ,.013,.008,-.070,.037],
        [(1,2),-.021,.014,.012,-.062,.039],
        [(1,3),-.023,.015,.013,-.066,.044],
        [(2,1),-.017,.013,.017,-.058,.048],
        [(2,2),-.025,.007,.010,-.069,.042],
        [(2,3),-.023,.006,.012,-.066,.040],
        [(3,1),-.012,.022,.023,-.055,.058],
        [(3,2),-.015,.019,.019,-.057,.052],
        [(3,3),-.020,.015,.015,-.061,.049]
        ]
        spec_table_aa={
        (1,1):[-.028,.012],
        (1,2):[-.025,.012],
        (1,3):[-.025,.013],
        (2,1):[-.018,.011],
        (2,2):[-.023,.005],
        (2,3):[-.023,.006],
        (3,1):[-.018,.017],
        (3,2):[-.020,.013],
        (3,3):[-.029,.011]}
        
        #table of volumes, row is site, with the first row being 1,
        #column is replicate, with leftmost being 1
        vol_table=[[128.,93.,93.],
                   [61.,44.,41.],  
                   [62.,52.,62.]]
        self.blanks=listToSpectrum([-.025,0,.004,-.07,.032]) 
        
        #now convert that table-like input array into spectrum objects
        self.sites={}#this will be a dictionary of replicates by site
        sitecount=len(spec_table)
        for i_line in range(sitecount):
            curline=spec_table[i_line]
            site,rep=curline[0]#get the site and replicate
            #turn all but the first column into a spectrum
            curspec=listToSpectrum(curline[1:])
            #get the after-acid data from a dictionary
            #keyed by (site,rep)
            curspec_aa=listToSpectrum(spec_table_aa[(site,rep)])
            #if a site doesn't exist, make a PigSpectSite object for it
            if site not in self.sites:
                self.sites[site]=PigSpectSite(site,self.blanks)

            #now make a replicate record for that site  
            self.sites[site].addRep(rep,{'ba':curspec,
                    'aa':curspec_aa,
                    'title':"Site %d, rep %d"%(site,rep),
                    'vol':vol_table[site-1][rep-1]})#use -1 to shift to 0-based
                        #array indices
            
            
    def __init__(self):
        """main routine for spectraphotometer CHL calculations"""
        self.sites={} #PigSpectSite items representing each site
        self.blanks={}#spectrophotometer blanks 
        self.getSpectrumInputs()#load the hardcoded input data
        for cursite in sorted(self.sites.keys()):
            self.sites[cursite].compute()#run the computations

    def describe(self):
        for cursite in sorted(self.sites.keys()):
            self.sites[cursite].describe()#show results on screen
            
    def makeAvgTable(self,out_fn):
        """Make an output table of the pigments averaged from the replicates"""
        heads="Site,Chl a,Chl a acid,Chl b,Chl c,Phaeopigment"
        form="%s,%.3f,%.3f,%.3f,%.3f,%.3f"
        tabout=TableMaker(out_fn,heads,form)
        for cursite in sorted(self.sites.keys()):
            site_obj=self.sites[cursite]
            ba=site_obj.avgpig_ba
            aa=site_obj.avgpig_aa
            tabout.writeRow([str(cursite),ba['chl_a'],aa['chl_a'],
                ba['chl_b'],ba['chl_c'],aa['phaeopig']])
        tabout.close()
        
#fluorometer section
TAU=2.27
#list of the valid values for sensitivity, with -1 being minimum sensitivity
VALID_SENS=[-1,3.16,10.0,31.6]
DOOR_FACTORS={
    False:{-1:0.08456, 3.16:0.02540, 10.0:0.00800, 31.6:0.00264},
    True:{ -1:0.000790,3.16:0.000252,10.0:0.000081,31.6:None}
}
def getDoorFactor(sens,is_x100):
    if sens not in VALID_SENS:#fuss if an invalid range setting is used
        raise ValueError("The input range setting is not one of the accepted values")
    dorf=DOOR_FACTORS[is_x100][sens]
    if dorf==None:
        raise ValueError("Door Factor not available for x100 x31.6")
    return dorf

def addToTup(intup,last_items):
    """Return a list containing whatever list or tuple is in 
    intup with the items in the list 'last_items' at the end."""
    return list(intup)+last_items
def make_FR_reps(reptups,acetone_vol,blanks):
    """Convenience function that makes a series of fluorRead objects 
    that share a standard acetone volume.  Reptups is a list of tuples
    where each tuple has the same variables in the same order as the 
    parameter list for the fluorRead constructor.  Acetone_vol 
    is the acetone volume shared between them.  The return value is a 
    list of fluorRead objects in the same order as the tuples that 
    defined them. """
    return [fluorRead(*addToTup(currargs,[acetone_vol,blanks])) for currargs in reptups]

class fluorRead:
    """This class represents a fluorometer reading"""
    def __init__(self,Fb,Fa,sens,is_x100,
                 filtvol,acetone_vol,blanks,title="Test"):
        self.title=title
        #take blanks into account, 
        blank_adj=blanks[is_x100][sens]
        self.Fb=Fb-blank_adj
        self.Fa=Fa-blank_adj

        self.sens=sens
        self.is_x100=is_x100 
        self.filtvol=filtvol
        self.acetone_vol=acetone_vol
        #self.chl is chl concentration
        #self.phaeo is phaeopigment concentration
        self.chl,self.phaeo=self.getChlAndPhaeo()
        #make the sensitivities into human readable strings
        if self.is_x100:
            xtag="x100"
        else:
            xtag="x1"
        if self.sens==-1:
            senstag='min sens'
        else:
            senstag='x%.2f'%self.sens
        #string describing self.sens, self.is_x100    
        self.sens_string="%s %s"%(xtag,senstag)
                
    def describe(self):
        print "  %s, Vol=%f"%(self.title,self.filtvol)
        print "  Fb=%f, Fa=%f"%(self.Fb,self.Fa)
        print "  sensitivity: %s"%(self.sens_string)
        print "  chl a: ", self.chl
        print "  phaeopigments: ",self.phaeo
        print 

    def getChlAndPhaeo(self):
        """Get the chl a and phaeopigment concentrations based 
        on before and after fluorescences(Fa,Fb), sensetivity settings and 
        volumes of acetone and sample that was filtered.  The
        return value is a tuple with the first item being chl_a and the second
        being phaeopigment concentration."""
        filtvol_L=self.filtvol*1e-3
        Fd=getDoorFactor(self.sens,self.is_x100)
        Fa=self.Fa;Fb=self.Fb#make shorthand to avoid clunky eq's
        #do shared parts of the equations
        taufac=TAU/(TAU-1)
        volfac=(self.acetone_vol/filtvol_L)
        
        chl_a=Fd*taufac*(Fb-Fa)*volfac
        phaeopig=Fd*taufac*(TAU*Fa-Fb)*volfac
        return (chl_a,phaeopig)
class PigFluoroSite:
    def __init__(self,sitenum):
        self.sitenum=sitenum
        self.reps={}
        
    def addRep(self,rep_index,inrep):
        """Add a replicate, inrep is a fluorRead object"""
        self.reps[rep_index]=inrep
    
    def makeAvg(self):
        """Make Chl and Phaeopigment averages, after all reps are in."""
        chl_list=[]
        phaeo_list=[]
        for i_currep in sorted(self.reps.keys()):
            phaeo_list.append(self.reps[i_currep].phaeo)
            chl_list.append(self.reps[i_currep].chl)
        self.chl_avg=np.mean(chl_list)
        self.phaeo_avg=np.mean(phaeo_list)        
    
    def describe(self):
        print "--------\nSite %d"%self.sitenum
        for i_currep in sorted(self.reps.keys()):
            self.reps[i_currep].describe()
        print "\n Averages: chl: %f, phaeo: %f"%(self.chl_avg,self.phaeo_avg)
        
class PigByFluoro:
    def __init__(self):
        self.sites={}
        #this dictionary represents a table where the 
        #key a tuple of (site,replicate), and the
        #value is a list representing columns of Fb, Fa, sens, and is_x100
        florins={
            (1,1):[4.0,2.1,3.16,False],
            (1,2):[4.0,2.1,3.16,False],
            (1,3):[3.9,2.1,3.16,False],
            (2,1):[2.2,1.4,10.0,False],
            (2,2):[2.1,1.5,10.0,False],
            (2,3):[2.2,1.5,10.0,False],
            (3,1):[5.9,3.2,10.0,False],
            (3,2):[6.0,3.4,10.0,False],
            (3,3):[6.0,3.6,10.0,False]
        }
        #this table of volumes is arrange with each row being a site,
        #and each column within being a replicate, the first of each is 1
        filt_vols=[[92.,88.,92.],
                   [46.,47.,52.],
                   [63.,60.,62.]]
        #acetone volumes grid, orginized like filt_vols
        ace_vols=[[6.,6.,6.],
                   [6.,6.1,6.],
                   [6.,6.,6.]]
        #define the blank measurements
        self.blanks={
            False:{-1:0.0,3.16:0.0,10.0:0.0,31.6:0.0},
            True:{-1:.1,3.16:.3,10.0:2.0,31.6:6.2}
        }
        #use the input tables to assemble a set of flourRead objects
        #organized by site and replicate
        for siterep in sorted(florins.keys()):
            meter_read=florins[siterep]
            sitenum,rep=siterep#split the site and replicate numbers
            #get the volumes for this replicate
            fvol=filt_vols[sitenum-1][rep-1]
            avol=ace_vols[sitenum-1][rep-1]
            #make a site for this site number if none exists
            if sitenum not in self.sites:
                self.sites[sitenum]=PigFluoroSite(sitenum)
            #make the replicate object
            nurep=fluorRead(meter_read[0],meter_read[1],
                meter_read[2],meter_read[3],fvol,avol,self.blanks,
                "Site %d Rep %d"%(sitenum,rep))
            self.sites[sitenum].addRep(rep,nurep)
        #compute all of the averages
        for sitenum,cursite in self.sites.items():
            cursite.makeAvg()
    def describe(self):
        for i_cursite in sorted(self.sites.keys()):
            self.sites[i_cursite].describe()
    def makeAvgTable(self,out_fn):
        """Make an output table of the pigments averaged from the replicates"""
        heads="Site,Chl a,Phaeopigment"
        form="%s,%.3f,%.3f"
        tabout=TableMaker(out_fn,heads,form)
        for cursite in sorted(self.sites.keys()):
            site_obj=self.sites[cursite]
            tabout.writeRow([str(cursite),site_obj.chl_avg,
                             site_obj.phaeo_avg])
        tabout.close()

def testFluorL5():
    blank_l5={
        False:{-1:0,3.16:0,10.0:0,31.6:0},
        True:{-1:.2,3.16:.6,10.0:1.9,31.6:5.8}
    }
    acetone_vol_l5=50.0
    #test_fluor=fluorRead(8.9,5.0,-1,True,250,acetone_vol)
    
    group_order=make_FR_reps([(8.9,5.0,-1,True,250)],acetone_vol_l5,blank_l5)
    #group_order=multiConstructFR([(8.9,5.0,-1,True,250,acetone_vol)])
    print "Lab 5 comparison test"
    for curgroup in group_order:
        curgroup.describe()

#main script
#run the spectral pigment analysis
spect_pigs=PigBySpectrum()
spect_pigs.makeAvgTable("chl_spect_avg.csv")
spect_pigs.describe()
fluoro_pigs=PigByFluoro()
fluoro_pigs.makeAvgTable("chl_fluoro_avg.csv")
fluoro_pigs.describe()
#test with the Lab5 flourometer data 
#testFluorL5()
