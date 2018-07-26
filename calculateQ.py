import csv
import matplotlib.pyplot as plt
import numpy 
import itertools
from scipy import stats

fileName = str(input())


global inp
inp = []
global out
out = []


def findMax(arr): ## Helper function for finding max frequency
    maxi = 0 
    index = 0
    for i in range(1,len(arr)/2):
        if arr[i] > maxi:
            index = i
            maxi = arr[i]
    return index

def findPreReqs(fileName,extent): ## Find pre-reqs (Sampling rate, Resonance Freq, Length of the file
    with open("details.csv") as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] == fileName:
                FirstElement = int(row[1])
                SamplingRate = float(row[2])
                Length = int(row[3])
                ResonanceFrequency = float(row[4])
                f.close()
                return [FirstElement,SamplingRate,Length,ResonanceFrequency]
    y = []
    count = 0
    with open(fileName) as f:
        FirstElement = 0
        reader = csv.reader(f)
        for row in reader:
            print row[0]
            if row[0][0] != '%':
                i = float(row[0])
                y.append(float(row[1]))
                count = count + 1
                break
            FirstElement = FirstElement + 1
        for row in itertools.islice(reader,0,1):
            j = float(row[0])
            y.append(float(row[1]))
            count = count + 1
        SamplingRate = 1/float((j-i))
        Length = 0
        for row in reader:
            if count < extent:
                y.append(float(row[1]))
                count = count+1
            Length = Length + 1
        f.close()
    Length = Length + FirstElement + 2
    y = numpy.abs(numpy.fft.fft(y))
    maxIndex = findMax(y)
    print maxIndex
    if maxIndex > len(y)/2:
        maxIndex = abs(maxIndex - len(y))
    ResonanceFrequency = SamplingRate*maxIndex/len(y)
    del y
    with open("details.csv",'ab') as f:
        writer = csv.writer(f)
        writer.writerow([fileName,FirstElement,SamplingRate,Length,ResonanceFrequency])
    return [FirstElement,SamplingRate,Length,ResonanceFrequency]

[FirstElement,SamplingRate,Length,ResonanceFrequency] = findPreReqs(fileName,10000)

def check(y,block,m):
    average = abs(sum(y)/block)
    if max(y)- average < 0.3*m: 
        return False
    return True

def humanIntervention(inp,out): ## Helper after calculating Q(Can select options to see what envelope is suitable)
    print "What do you want?"
    x = []
    y = []
    while 1:
        z = str(input("Next"))
        x = []
        y = []
        if z == "Plot":
            plt.plot(inp,out)
            plt.show()
            print len(inp)
        elif z == "Envelope":
            envelope(out,y,inp,x)
            out = y[:]
            inp = x[:]
            y = []
            x = []
        elif z == "Abs":
            out = numpy.fft.fft(out)
            out[0] = 0
            out = list(numpy.abs(numpy.fft.ifft(out)))
        elif z == "Q":
            p = numpy.polyfit(inp,numpy.log(out),1)
            plt.plot(inp,numpy.log(out),'r')
            y = numpy.multiply(inp,p[0])
            y = numpy.add(y,p[1])
            plt.plot(inp,y,'b')
            plt.show()
            k = str(input())
            if k == "Fine":
                return float(-1*ResonanceFrequency*22)/(7*p[0])
        else:
            continue
                
def calculateQ(inp,out,block): ## Calculates Q. Block is the number of data items taken one at a time.
    with open(fileName) as f:
        reader = csv.reader(f)
        start = FirstElement
        x = []
        y = []
        for row in itertools.islice(reader,start,start+block):
            x.append(float(row[0]))
            y.append(float(row[1]))
        if(len(y)!=0):
            fft = numpy.abs(numpy.fft.fft(y))
            maxIndex = findMax(fft)
            maxIndex = int(len(y)/maxIndex)
            i = 0
            while i < len(y):
                inp.append(x[i])
                out.append(y[i])
                i = i + maxIndex
        count = start+block
        m = y[0] - abs(sum(y)/block)
        while count < Length:
            x = []
            y = []
            for row in itertools.islice(reader,0,start+block):
                x.append(float(row[0]))
                y.append(float(row[1]))
            if(len(y)!=0):
                if check(y,block,m) is False:
                    break
                fft = numpy.abs(numpy.fft.fft(y))
                maxIndex = findMax(fft)
                maxIndex = int(len(y)/maxIndex)
                i = 0
                while i < len(y):
                    inp.append(x[i])
                    out.append(y[i])
                    i = i + maxIndex
            count = count + block
    out = numpy.fft.fft(out)
    out[0] = 0
    out = numpy.abs(numpy.fft.ifft(out))
    return humanIntervention(inp,out)

def envelope(data,derivative,inp,values):  ## Takes envelope by selecting points where derivative changes sign. This is used in the calculateQ
    derivative.append(data[0])
    values.append(inp[0])
    for i in range(1,len(data)-1):
        if data[i]-data[i-1]>=0 and data[i]-data[i+1]>=0:
            derivative.append(data[i])
            values.append(inp[i])
    derivative.append(data[len(data)-1])
    values.append(inp[len(data)-1])
    return    
    
def calculateCorrelation(fileName,extent): ## Calculate correlation between x and y
    with open(fileName) as f:
        reader = csv.reader(f)
        start = FirstElement
        inp = []
        out = []
        for row in itertools.islice(reader,start,start+extent):
            inp.append(float(row[1]))
            out.append(float(row[2]))
        return stats.linregress(inp,out)

def analyseY(fileName,extent): ## Extract the y data to a certain extent
    inp = []
    out = []
    with open(fileName) as f:
        reader = csv.reader(f)
        start = FirstElement
        for row in itertools.islice(reader,start,start+extent):
            inp.append(float(row[0]))
            out.append(float(row[2]))
        return [inp,out]

def plotFrequencies(fileName,block): ## Plots the frequencies in x and y data for blocks against each other
    x_freqs = []
    y_freqs = []
    with open(fileName) as f:
        reader = csv.reader(f)
        start = FirstElement
        x = []
        y = []
        for row in itertools.islice(reader,start,start+block):
            x.append(float(row[1]))
            y.append(float(row[2]))
        if(len(x)!=0):
            fft = numpy.abs(numpy.fft.fft(x))
            maxIndex = findMax(fft)
            x_freqs.append(maxIndex*SamplingRate/block)
            fft = numpy.abs(numpy.fft.fft(y))
            maxIndex = findMax(fft)
            y_freqs.append(maxIndex*SamplingRate/block)
        count = start+block
        while count < Length:
            x = []
            y = []
            for row in itertools.islice(reader,start,start+block):
                x.append(float(row[1]))
                y.append(float(row[2]))
            if(len(x)!=0):
                fft = numpy.abs(numpy.fft.fft(x))
                maxIndex = findMax(fft)
                x_freqs.append(maxIndex*SamplingRate/block)
                fft = numpy.abs(numpy.fft.fft(y))
                maxIndex = findMax(fft)
                y_freqs.append(maxIndex*SamplingRate/block)
            count = count + block
    plt.figure(1)
    plt.plot(x_freqs)
    plt.figure(2)
    plt.plot(y_freqs)
    plt.figure(3)
    plt.plot(x_freqs,y_freqs)
    plt.show()
    print len(x_freqs)
    print len(y_freqs)
    
    

    
