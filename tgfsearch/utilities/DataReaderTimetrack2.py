#Version 230226
#Changed sign of unixtimecorrection to subtraction in processdatatiming()
#Changed location of best time tag string in thor data to lines[i-4] (first tag) in thorFileToData()
#Added ucsc_timestring_to_datetime() routine to correctly pad microseconds in UCSC time tags with 0s on the left,
#         needed by datetime.strptime()
#Added multifilesToData() routine to read multiple files into a single dataset
#Added warnings about unexpected behavior across frame boundaries to thorFileToData()
#5/6/24: added handling of the case where there are columns for specific flags (like 'pps') instead of a
#        'flags' column like in normal THOR data in getdatafromlmthor()
#5/6/24: Changed "if mode is" to "if mode ==" on advice of new warnings.
#7/5/24 Fixed condition where there is a rollover between the last PPS of the prior frame and the end
#       of the of the prior frame, which was messing up all the counts up to the first pps of the new frame.

import re # Regular expressions module for string operations
import pdb  #syntax to use is pdb.set_trace() for equivalent to an IDL "stop"
import json
import pandas as pd
import glob
import os
import numpy as np
import scipy
from scipy.stats import poisson
from datetime import datetime
import gzip
from matplotlib import pyplot as plt

def multifilesToData(filefile):
    f=open(filefile,'r')
    lines = f.read().splitlines()
    started=0
    passtime = {"lastsod":-1.0, "ppssod":-1.0, "lastunix":-1.0, "ppsunix":-1.0, "lastwc":-1, "ppswc":-1, "hz":8e7, "started":0}
    for s in lines:
        print("using file: "+s)
        data,passtime = fileNameToData(s,passtime)
        if (started==0):
            alldata = data
            started=1
        else:
            #pdb.set_trace()
            alldata=pd.concat([alldata,data],axis=0)
    return(alldata)

def ucsc_timestring_to_datetime(timeString):
    parts = timeString.split()
    microsec = parts[-1]
    microsec = microsec.rjust(6,'0')
    newString = parts[0]+' '+parts[1]+' ' +parts[2]+' '+parts[3]+' '+parts[4]+' '+parts[5]+' '+microsec
    headerDT = datetime.strptime(newString, "%Y %m %d %H %M %S %f")
    return(headerDT)

def fileNameToData(fileName,passtime):
    if(fileName[-2:] == 'gz'):        # This decompresses the .gz file if it's gz, that way you don't have to uncompress a file first if you don't want to. Works either way.
        f = gzip.open(fileName, 'rt')
    else:
        f = open(fileName)
    lines = f.read().splitlines()
    if ("xtr" in fileName):         # This means it's a trace file, so go to the trace reading function
        return(traceFileToData(lines))
    if(len(lines[0]) == 3):         # Only the thor LM files start with a 3 character (including a space) line, so that's the fastest way to identify that format
        data,passtime = thorFileToData(lines,passtime)
    else:
         mode = 1
         if len(lines[2]) < 50: #Early LM files (mode 1) Just alternate (Time) - (Buffer) - (Time) - (Buffer) etc. 
                                #The 2nd version has 4 time tags, so seeing that the second line is a time tag and not a full buffer tells you the type.
            mode = 2
         data = lmFileToData(lines, mode)
    #print('! ',fileName)
    #print(data)
    #print(passtime)
    return(data,passtime)


##Perhaps these next two files could be combined, they're very similar, the main difference is the call to "getDataFromLMThor" and the different numbers used for iteration
def thorFileToData(lines, passtime):
    dataList = list()
    for i in range(5, len(lines), 6):
        lmString = lines[i]
        timeString = lines[i - 4] #John had this wrong (he had i-2, which = 3rd stamp)
        headerDT = ucsc_timestring_to_datetime(timeString)

        data = getDataFromLMthor(lmString)
        data, newpasstime = processDataTiming(data, passtime, headerDT, mode=0)

        #In this condition, you are near the end of the day and the buffer goes into the next day.
        w = np.where( np.array(data['SecondsOfDay']) < 0.) 
        mn = max( np.array(data['SecondsOfDay']))
        if (len(w[0]) > 0) and (mn < 5000.):
            data['SecondsOfDay'] = data['SecondsOfDay'] + 86400.
        #What happens if we just keep going, though?

        tfirst = data['SecondsOfDay'][0]
        tlast =  data['SecondsOfDay'][len(data)-1]

        if (i > 5):
            dt=tfirst-passtime['lastsod']
            if (dt > 0.5):
                print('Long gap (>0.5s) between frames at: ',tfirst)
                #pdb.set_trace()
            if (dt < 0.0):
                print('Anomalous clock backwards at: ',tfirst,passtime['lastsod'], dt)
#                pdb.set_trace()
#            if (dt < -50000.):
#                print('...due to end of day.  Fixing.')
                #pdb.set_trace()
#                w = np.where( np.array(data['SecondsOfDay']) - tfirst < -50000.) 
 #               pdb.set_trace()
                #data['SecondsOfDay'][w] =  data['SecondsOfDay'][w] + 86400.000

        passtime=newpasstime
        dataList.append(data)
    data = pd.concat(dataList)
    return(data,passtime)


def lmFileToData(lines, mode):
    dataList = list()
    if mode== 1:
        start = 2
        increment = 2
        tStamp = 1
    elif mode== 2:
        start = 7
        increment = 6
        tStamp = 0 #NEED TO DOUBLE CHECK THIS WITH JEFF (was 2 -- experimenting)
    for i in range(start, len(lines), increment):
        lmString = lines[i]
        timeString = lines[i + tStamp]
        headerDT = ucsc_timestring_to_datetime(timeString)
        
        data = getDataFromLM(lmString, mode)
        if (i + increment < len(lines)):
            nextLmString = lines[i + increment]
            nextData = getDataFromLM(nextLmString, mode)
            match = data['wc',data.index[0]] == nextData['wc',data.index[0]]
            if (match):
                print("duplicate buffer " + str(i))
                continue
        data,passtime = processDataTiming(data, passtime,headerDT, mode=mode)
        dataList.append(data)
    data = pd.concat(dataList)
    return(data)


def getDataFromLMthor(lmString):
    jsonDict = json.loads(re.sub("eRC[0-9]{4} ", "", lmString))['lm_data']
    data = pd.DataFrame.from_dict(jsonDict)
    if 'pps' in data.columns:
        data['PPS'] = data['pps'].copy()
    else:
        data['PPS'] = pd.Series([int('{0:08b}'.format(x)[-8]) for x in data['flags']]) ## very useful column for future operations - separates out GPS signal from the Flags column
    data['gpsSync'] = False
    data['UnixTime']=0.
    data['SecondsOfDay']=0.

    return(data)


def processDataTiming(data, passtime, headerDT, mode):

    lst = len(data)-1
    lstpps = -1

    #Get rid of any rollovers within the frame

    rolloverLength = pow(2, 36)
    if mode == 1:
        rolloverLength = pow(2, 32)
    elif mode == 2:
        rolloverLength = pow(2, 48)
    rolloverCorrection = (data['wc'].diff() < -rolloverLength/4.).cumsum() * rolloverLength  
    ## Every time you encounter the wallclock going significantly backwards (1/4 of range), 
    ## assume it's rollover. Count how many times its rolled over and multiply that by 2^36 for each row's correction.

    data['wc'] = data['wc'] + rolloverCorrection ## Changing the contents of the object passed in the function - 
                                                 ##there might be a better way to do this but it should work since 
                                                 ##Python passes by object reference.

    #But what if there was a rollover between the last PPS of the prior frame and the start of this frame???? -- Note to add code for that. :p
    #7/5/24 Here is is, looks like it works (?)
    lastpps_wc = passtime['ppswc']
    first_wc_here = data['wc'][0]
    if (first_wc_here - lastpps_wc < -rolloverLength/4.):
        passtime['ppswc'] -= rolloverLength

    
    #Get the header timestamp in unix time, and unix time at start of second, and start of day
    header_unix = headerDT.timestamp()
    header_unix_trunc = int(header_unix)*1.0
    daystart_unix = datetime(headerDT.year, headerDT.month, headerDT.day).timestamp() ## subtracts out 00:00:00 of the earliest day in UNIX time.         

    #Get PPS events (if available)

    if (mode==0):
        pps = np.where(data['PPS'])
        pps = pps[0]
        npps = len(pps)
    else:
        npps = 0

    #Create a new time structure.  It will get updated as needed and then pushed back into "passtime"
    newtime=passtime.copy()

    #First method: in-frame only, last count = timetag:
    m1_unix =  header_unix - (data['wc'][lst] - data['wc'])/newtime['hz']

    #Second method: across-frame, not using PPS:
    if ( passtime['started']>0 ):
        m2_unix =  passtime['lastunix']  + (data['wc'] - newtime['lastwc']) / newtime['hz']
    else:
        m2_unix = m1_unix*0. - 1.

    #If it's the very first buffer, base everything on the header for now, but it will get fixed later
    #Last event in frame is assumed to be exactly the header time.  Wall clock rate is assumed to be precise (hz=8e7, the default)
    if (npps > 0):
        lstpps = pps[-1]
            
    if npps < 2:
        if ( passtime['started']==0 ):
            if (mode==0):
                print('THOR data begins with a frame with 0 or 1 PPS!')
            data['UnixTime'] = m1_unix
        else:
            data['UnixTime'] = m2_unix
            print('THOR data has a frame with ' ,npps,  ' PPS: ',data['wc'][0] )

    else:  #This is THOR data with PPS flags in it

        #Find the frequencies in between each pair of PPS:
        frequencies = m1_unix*0. + 8e7
#        if ( passtime['started']>0 ):
#              pdb.set_trace()
        for i in np.arange(len(pps)-1):
            testhz = data['wc'][pps[i+1]] - data['wc'][pps[i]]

            if abs(testhz - 8e7) < 1e4:
                frequencies[pps[i]:pps[i+1]+1] = testhz
                data.loc[   pps[i]:pps[i+1]+1,'gpsSync'] = True
            else:  #if not synced, try guessing that one PPS was skipped; otherwise out of sync.
                if (abs(testhz/2.0 - 8e7) < 1e4):
                    frequencies[pps[i]:pps[i+1]+1] = testhz/2.0
                    data.loc[pps[i]:pps[i+1]+1,'gpsSync'] = True
                    print('Missing one PPS: ',headerDT, testhz, testhz/2.0)
                else:
                    data.loc[pps[i]:pps[i+1]+1,'gpsSync'] = False
                    print('GPS sync failed.  Factor: ', testhz/8.0e7, ' at: ',headerDT)
            if i==0:
                frequencies[0:pps[0]] = testhz
            if i==len(pps)-2:
                frequencies[pps[i+1]:] = testhz

        #If prior PPS is available from before the frame, improve the starting events before first PPS:
        if ( passtime['started']>0 ):
            testhz2 = data['wc'][pps[0]] - passtime['ppswc']
            if abs(testhz2 - 8e7) < 1e4:
                frequencies[0:pps[0]] = testhz2
                data.loc[0:pps[0],'gpsSync'] = True
            else:  #if not synced, try guessing that one PPS was skipped; otherwise out of sync.
                if (abs(testhz2/2.0 - 8e7) < 1e4):
                    frequencies[0:pps[0]] = testhz2/2.0
                    data.loc[0:pps[0] ,'gpsSync'] = True
                    print('Missing one PPS across frame boundary: ',headerDT,testhz2,testhz2/2.0)
                else:
                    data.loc[0:pps[0],'gpsSync'] = False
                    print('GPS sync failed.  Factor: ', testhz2/8.0e7, ' at: ',headerDT)
            
        #Extrapolate the events after the last PPS:
        frequencies[pps[-1]:lst] = frequencies[pps[-1]-1]

        #Now get full array of unix time by two methods.
        #Start with the reverse method (the only one when it's the very first buffer):

        index = 0               #index==0 is the *last* pps 
        m3_unix = m1_unix*0.
        rev_pps = np.flip(pps)  #rev_pps[index] = i   as used below

        for i in rev_pps:

            if (index == 0):             #Last pps in frame, first to be examined
                m3_unix[i] = int(m1_unix[i]+0.1)
                m3_unix[i:lst+1] = m3_unix[i] + (data['wc'][i:lst+1] - data['wc'][i])/frequencies[i:lst+1]
                testdiff = m3_unix[i] - m1_unix[i]
            else:                        #all but the very first pps
                known = rev_pps[index-1]
                m3_unix[i:known] = m3_unix[known] - (data['wc'][known] - data['wc'][i:known])/frequencies[i:known]
            if (index==npps-1):                        #first pps in frame, last to be examined
                m3_unix[0:i] =  m3_unix[i] - (data['wc'][i] - data['wc'][0:i])/frequencies[0:i]
            index += 1

        #Forward method if available:
        #Find the PPS events & walk through each interval, starting with the one that connects to the previous 
        #buffer through "passtime":

        m4_unix = m1_unix*0. - 1.

        if (passtime['started']>0 ):

            index = 0

            for i in pps:

                m4_unix[i] = int(m2_unix[i])  

                if (index == 0):             #first pps in frame, first to be examined
                    m4_unix[0:i+1] =  passtime['ppsunix'] + (data['wc'][0:i+1] - passtime['ppswc'])/frequencies[0:i+1]

#                   m4_unix[0:i] = m4_unix[i] - (data['wc'][i] - data['wc'][0:i])/frequencies[0:i]

                if (index < len(pps)-1):     #all but the very last pps
                    m4_unix[i:pps[index+1]] = m4_unix[i] + ( data['wc'][i:pps[index+1]] - data['wc'][i])/frequencies[i:pps[index+1]]
                else:                        #last pps in frame, last to be examined
                    m4_unix[i:lst+1] = m4_unix[i] + (data['wc'][i:lst+1] - data['wc'][i])/frequencies[i:lst+1]

                index = index + 1

        if ( passtime['started']==0 ):
            data['UnixTime'] = m3_unix
        else:
            diff = m3_unix-m4_unix
            #if ( (min(diff) < -0.5) or (max(diff) > 0.5) ):
            #    pdb.set_trace()
            data['UnixTime'] = m4_unix


####################################################################################

        #Check last pps against integer of header second:
        #diff = header_unix_trunc - newtime['ppsunix']
        #if (abs(diff) > 1e-2 and passtime['started'] > 0):
        #    print("Warning: mismatch for integer second of last PPS. Header: ",header_unix,' PPS time: ',newtime['ppsunix'])
        #    if (abs(header_unix - header_unix_trunc) > 0.1):
        #        print("    WARNING: header time is NOT near a second boundary.  Using extrapolated solution.")
        #        #pdb.set_trace()

    #Every time after the very first time this routine is called, the "started" flag will = 1.
    data['SecondsOfDay'] = data['UnixTime'] - daystart_unix
    newtime['lastunix'] = data['UnixTime'][lst]
    newtime['lastsod'] = data['SecondsOfDay'][lst]
    newtime['lastwc'] = data['wc'][lst] % rolloverLength
    if (npps > 0): 
        newtime['ppswc'] = data['wc'][pps[-1]] % rolloverLength
        newtime['ppsunix'] = data['UnixTime'][lstpps]
        newtime['ppssod'] = data['SecondsOfDay'][lstpps]

    #if (newtime['started']==1):
    #    f = open("test.txt", "w")
    #    for i in np.arange(lst+1):
    #        f.write("%s %s %s %s %s %s\n" %(i, data['wc'][i], m1_unix[i] - daystart_unix, m2_unix[i] - daystart_unix, m3_unix[i] - daystart_unix, m4_unix[i] - daystart_unix))
    #    f.close()
    #    pdb.set_trace()

    newtime['started'] = 1


    passtime = newtime
 
    return(data,passtime)


def getDataFromLM(lmString, mode):
    if mode== 1:
        start = 2
        rowLength = 3
    elif mode== 2:
        start = 9
        rowLength = 6
    lmStrings = lmString.split(" ")[start:]
    buffer = [int(y) for y in lmStrings]
    data = np.reshape(buffer, (int(len(buffer)/rowLength), rowLength))
    data = pd.DataFrame(data)
    #data = pd.concat(data, ignore_index = True)
    data = data[~data.duplicated()]
    coarsetick = 65536
    if mode== 2:
        wc = (data[2] + data[3] * coarsetick + data[4] * coarsetick * coarsetick) 
        energy = data[1]
        ticks = pd.Series([int('{0:08b}'.format(x)[-8]) for x in data[5]])
        flags = data[5]
    elif mode== 1:
        wc = (data[1] + data[2] * coarsetick) 
        energy = data[0]
        ticks = data[0] * 0 ## For some reason makes the pd.concat phase much faster than doing (repeat(0, len(data)))
        flags = data[0] * 0
    newData = pd.concat([energy.rename('energy'), wc.rename('wc'), ticks.rename('PPS'), flags.rename('flags')], axis = 1)
    return(newData)


def traceFileToData(lines):
    dataList = list()
    for i in range(4, len(lines), 5):
            dataLine = lines[i]
            timeString = lines[i - 1]
            jsonDict = json.loads(re.sub("eRC[0-9]{4} [0-9]", "", lines[i]))
            data = pd.DataFrame.from_dict(jsonDict)
            data['BufferNo'] = int(lines[i][8:9]) # There should be some significance to how this affects the time relationships?
            data['DateTime'] = datetime.strptime(timeString, "%Y %m %d %H %M %S %f")
            data['Seconds'] = [x * 1.25e-8 for x in range(len(data))]
            dataList.append(data)
    newData = pd.concat(dataList)
    return(newData)




    
