#Version 230226
#Changed sign of unixtimecorrection to subtraction in processdatatiming()
#Changed location of best time tag string in thor data to lines[i-4] (first tag) in thorFileToData()
#Added ucsc_timestring_to_datetime() routine to correctly pad microseconds in UCSC time tags with 0s on the left,
#         needed by datetime.strptime()
#Added multifilesToData() routine to read multiple files into a single dataset
#Added warnings about unexpected behavior across frame boundaries to thorFileToData()

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
    passtime = {"lastsod":-1.0, "ppssod":-1.0, "lastunix":-1.0, "ppsunix":-1.0, "lastwc":0, "ppswc":0, "hz":8e7, "started":0}
    for s in lines:
        #print("using file: "+s)
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
    return(data,passtime)


##Perhaps these next two files could be combined, they're very similar, the main difference is the call to "getDataFromLMThor" and the different numbers used for iteration
def thorFileToData(lines, passtime):
    dataList = list()
    for i in range(5, len(lines), 6):
        lmString = lines[i]
        timeString = lines[i - 4] #John had this wrong (he had i-2, which = 3rd stamp)
        headerDT = ucsc_timestring_to_datetime(timeString)

        data = getDataFromLMthor(lmString)
        data, passtime = processDataTiming(data, passtime, headerDT, mode=0)
        tfirst = data['SecondsOfDay'][0]
        tlast =  data['SecondsOfDay'][len(data)-1]
 
        if (i > 5):
            dt=tfirst-tlast_prev
            if (dt > 0.5):
                #print('Long gap (>0.5s) between frames at: ',tfirst)
                pass
            if (dt < 0.0):
                #print('Anomalous clock backwards at: ',tfirst,tlast_prev, dt)
                pass

        tlast_prev = tlast
        dataList.append(data)
    data = pd.concat(dataList)
    return(data,passtime)


def lmFileToData(lines, mode):
    dataList = list()
    if mode is 1:
        start = 2
        increment = 2
        tStamp = 1
    elif mode is 2:
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
                #print("duplicate buffer " + str(i))
                continue
        data,passtime = processDataTiming(data, passtime,headerDT, mode=mode)
        dataList.append(data)
    data = pd.concat(dataList)
    return(data)


def getDataFromLMthor(lmString):
    jsonDict = json.loads(re.sub("eRC[0-9]{4} ", "", lmString))['lm_data']
    data = pd.DataFrame.from_dict(jsonDict)
    data['PPS'] = pd.Series([int('{0:08b}'.format(x)[-8]) for x in data['flags']]) ## very useful column for future operations - separates out GPS signal from the Flags column
    data['gpsSync'] = False
    data['UnixTime']=0.
    data['SecondsOfDay']=0.

    return(data)


def processDataTiming(data, passtime, headerDT, mode):

    lst = len(data)-1

    #Get rid of any rollovers within the frame

    rolloverLength = pow(2, 36)
    if mode == 1:
        rolloverLength = pow(2, 32)
    elif mode == 2:
        rolloverLength = pow(2, 48)
    rolloverCorrection = (data['wc'].diff() < -rolloverLength/4.).cumsum() * rolloverLength  ## Every time you encounter the wallclock going significantly backwards (1/4 of range), 
                                                                        ## assume it's rollover. Count how many times its rolled over and multiply that by 2^36 for each row's correction.
    data['wc'] = data['wc'] + rolloverCorrection ## Changing the contents of the object passed in the function - there might be a better way to do this but it should work since 
                                                 ##Python passes by object reference.

    #But what if there was a rollover between the last PPS of the prior frame and the start of this frame???? -- Note to add code for that. :p
    

    #Get the header timestamp in unix time, and unix time at start of second, and start of day
    header_unix = headerDT.timestamp()
    header_unix_trunc = int(header_unix)*1.0
    daystart_unix = datetime(headerDT.year, headerDT.month, headerDT.day).timestamp() ## subtracts out 00:00:00 of the earliest day in UNIX time.         

    #Get number of PPS events (if available)
    if (mode==0):
        npps = np.sum(data['PPS'])
    else:
        npps = 0

    #Create a new time structure.  It will get updated as needed and then pushed back into "passtime"
    newtime=passtime.copy()

    #When no PPS and it's the very first buffer, use the header time to create timing for the whole buffer.
    #Last event in frame is assumed to be exactly the header time.  Wall clock rate is assumed to be precise (hz=8e7, the default)
    if ( (newtime['started']==0 and npps==0 and mode==0) or (newtime['started']==0 and mode > 0) ):
        if (mode==0):
            #print('THOR data begins with a frame with no PPS!')
            pass
        data['UnixTime'] = header_unix - (data['wc'][lst] - data['wc'])/newtime['hz']
        data['SecondsOfDay'] = data['UnixTime'] - daystart_unix
        newtime['lastunix'] = data['UnixTime'][lst]
        newtime['lastsod'] = data['SecondsOfDay'][lst]
        newtime['lastwc'] = data['wc'][lst]

    #If it's not the first buffer but there are no PPS (whether because there can't be in this mode, or whether missing),
    #extend the prior timing solution through the frame:
    elif npps==0:

        timediffs = (data['wc'] - newtime['lastwc']) / newtime['hz']
        data['SecondsOfDay'] = newtime['lastsod'] + timediffs
        data['UnixTime'] = newtime['lastunix']  + timediffs
        newtime['lastunix'] = data['UnixTime'][lst]
        newtime['lastsod'] = data['SecondsOfDay'][lst]
        newtime['lastwc'] = data['wc'][lst]

    else:

        #Find the PPS events & walk through each interval, starting with the one that connects to the previous 
        #buffer through "passtime":
        pps = np.where(data['PPS'])
        pps = pps[0]
        index = 0

        for i in pps:
            if (index==0):
                start = 0
            else:
                start = pps[index-1]+1

            #Update the true wall clock frequency if possible, and flag whether it's synced to GPS:
            if (passtime['started'] > 0 or index > 0):
                testhz = data['wc'][i] - newtime['ppswc']
                if abs(testhz - 8e7) < 1e4:
                    newtime['hz'] = testhz
                    data.loc[start:i+1,'gpsSync'] = True
                else:  #if not synced, try guessing that one PPS was skipped; otherwise out of sync.
                    if (abs(testhz/2.0 - 8e7) < 1e4):
                        newtime['hz'] = testhz/2.0
                        data.loc[start:i+1,'gpsSync'] = True
                        #print('Missing one PPS: ',headerDT)
                    else:
                        data.loc[start:i+1,'gpsSync'] = False
                        #print('GPS sync failed.  Factor: ', testhz/8.0e7, ' at: ',headerDT)
                    

            #Calculate the Unix time and seconds of day up through this PPS count.
            #"last" means either last event of the prior frame, if this is the stretch
            #before the first PPS of the buffer; otherwise it means the prior PPS event.
            timediffs = (data['wc'][start:i+1] - newtime['lastwc']) / newtime['hz']
            data.loc[start:i+1,'SecondsOfDay'] = newtime['lastsod'] + timediffs
            data.loc[start:i+1,'UnixTime'] = newtime['lastunix']  + timediffs

            newtime['ppssod'] = data['SecondsOfDay'][i]
            newtime['ppsunix'] = data['UnixTime'][i]
            newtime['ppswc'] = data['wc'][i]

            index = index+1

        #Fill in the last counts after the last PPS:
        timediffs = (data['wc'][i+1:] - newtime['ppswc']) / newtime['hz']
        data.loc[i+1:,'SecondsOfDay'] = newtime['ppssod'] + timediffs
        data.loc[i+1:,'UnixTime'] = newtime['ppsunix']  + timediffs       

        #If this is the very first buffer, lock the last PPS to the truncated header time:
        if (newtime['started']==0):
            fixtime = header_unix_trunc - data['UnixTime'][i]
            data['UnixTime'] = data['UnixTime'] + fixtime
            data['SecondsOfDay'] = data['UnixTime'] - daystart_unix

        #Set last event parameters for passtime (ppswc etc.. are already set)
        newtime['lastsod'] = data['SecondsOfDay'][lst]
        newtime['lastunix'] = data['UnixTime'][lst]
        newtime['lastwc'] = data['wc'][lst]

        #Check last pps against integer of header second:
        diff = header_unix_trunc - newtime['ppsunix']
        if (abs(diff) > 1e-2 and passtime['started'] > 0):
            #print("Warning: mismatch for integer second of last PPS. Header: ",header_unix,' PPS time: ',newtime['ppsunix'])
            if (abs(header_unix - header_unix_trunc) > 0.1):
                #print("    WARNING: header time is NOT near a second boundary.  Using extrapolated solution.")
                #pdb.set_trace()
                pass

    #Every time after the very first time this routine is called, the "started" flag will = 1.
    passtime = newtime
    passtime['started']=1
    passtime['lastwc'] = passtime['lastwc'] % rolloverLength
    passtime['ppswc'] = passtime['ppswc'] % rolloverLength

    return(data,passtime)


def getDataFromLM(lmString, mode):
    if mode is 1:
        start = 2
        rowLength = 3
    elif mode is 2:
        start = 9
        rowLength = 6
    lmStrings = lmString.split(" ")[start:]
    buffer = [int(y) for y in lmStrings]
    data = np.reshape(buffer, (int(len(buffer)/rowLength), rowLength))
    data = pd.DataFrame(data)
    #data = pd.concat(data, ignore_index = True)
    data = data[~data.duplicated()]
    coarsetick = 65536
    if mode is 2:
        wc = (data[2] + data[3] * coarsetick + data[4] * coarsetick * coarsetick) 
        energy = data[1]
        ticks = pd.Series([int('{0:08b}'.format(x)[-8]) for x in data[5]])
        flags = data[5]
    elif mode is 1:
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




    
