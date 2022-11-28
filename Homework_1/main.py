import numpy as np
import scipy.io.wavfile as wf
from tkinter import *
from tkinter.ttk import *
from tkinter import messagebox
import matplotlib.pyplot as plt



class App(Frame):
    def __init__(self, master=None):
        Frame.__init__(self, master)
        Style().theme_use('default')
        self.master = master
        self.master.title("Frequency Spectrum Analyzer")

        self.flag = 0

        self.path = StringVar(self, value="data/ribizlaM.wav")
        self.path2 = StringVar(self, value="data/slikaW.wav")
        Label(self, text="Path of first fIle").grid(row=0, column=0, padx=10, pady=10)
        Entry(self, textvariable=self.path).grid(row=0, column=1, padx=10, pady=10)

        Label(self, text="Path of second fIle").grid(row=1, column=0, padx=10, pady=10)
        Entry(self, textvariable=self.path2).grid(row=1, column=1, padx=10, pady=10)

        self.win_size = IntVar(self, value=10)
        Label(self, text="Window size (ms)").grid(row=2, column=0, padx=10, pady=10)
        Entry(self, textvariable=self.win_size).grid(row=2, column=1, padx=10, pady=10)



        self.win_funs = {"Hanning":0, "Hamming":1, "None":2}
        self.win_fun = IntVar(self, value=0)
        Label(self, text="Windowing").grid(row=3, column=0, padx=10, pady=10)
        for i, (k, v) in enumerate(self.win_funs.items()):
            Radiobutton(self, variable=self.win_fun, text=k, value=v).grid(row=4, column=i, padx=10, pady=10)

        Button(self, command=self.openFile, text="Load WAV").grid(row=5, column=0, padx=10, pady=10)
        Button(self, command=self.cutFile, text="Cut First WAV").grid(row=5, column=1, padx=10, pady=10)
        Button(self, command=self.drawBoth, text="Cut Both WAV").grid(row=5, column=2, padx=10, pady=10)
        Button(self, command=self.histogram, text="Show histogram").grid(row=6, column=0, padx=10, pady=10)
        Button(self, command=self.spectogram, text="Show spectrogram").grid(row=6, column=1, padx=10, pady=10)

        self.pack(fill=BOTH, expand=1)

    def openFile(self):

        try:
            self.samplerate, self.data = wf.read(self.path.get())
            self.samplerate2, self.data2 = wf.read(self.path2.get())

        except:
            messagebox.showinfo("Error", "File not found")


    def cutFile(self):



        #print(self.samplerate)

        self.T = [i / self.samplerate for i in range(0, len(self.data))]#vreme izmedju svakog pojedinacnog sempla
        self.T2 = [i / self.samplerate2 for i in range(0, len(self.data2))]

        #print(self.T)

        noiseArea = int(self.samplerate * 0.1)
        noise = np.abs(self.data[:noiseArea])
        self.L = np.mean(noise) + 2 * np.std(noise)


        noiseArea2 = int(self.samplerate2 * 0.1)
        noise2 = np.abs(self.data2[:noiseArea2])
        self.L2 = np.mean(noise2) + 2 * np.std(noise2)

        # L - average noise
        # N - number of samples in window
        self.N = int(self.samplerate * 0.01)
        print(self.N)

        self.length = len(self.data)

        self.words = [1 if np.mean(np.abs(self.data[i:i+self.N])) > self.L else 0 for i in range(0, self.length, self.N)]



        self.N2 = int(self.samplerate2 * 0.01)

        self.length2 = len(self.data2)

        self.words2 = [1 if np.mean(np.abs(self.data2[i:i + self.N2])) > self.L2 else 0 for i in range(0, self.length2, self.N2)]
        


        self.flattenup(self.words, 12)

        self.flattendown(self.words, 12)

        
        
        self.flattenup(self.words2, 12)

        self.flattendown(self.words2, 12)



        try:



            self.rec = self.drawBorders(self.words, self.data, self.T, self.N)

        except:
            messagebox.showinfo("Error", "Noise detected")
            self.drawNoise()


        self.window_function()

    def drawBoth(self):
        start1 = self.words.index(1)
        end1 = len(self.words) - 1 - self.words[::-1].index(1) + 1

        start2 = self.words2.index(1)
        end2 = len(self.words2) - 1 - self.words2[::-1].index(1) + 1

        f = plt.figure()
        plt.plot(self.T[::], self.data[::], 'b')
        plt.axvline(x=self.T[start1 * self.N], color='r')
        plt.axvline(x=self.T[end1 * self.N], color='r')

        plt.plot(self.T2[::], self.data2[::], 'y')
        plt.axvline(x=self.T2[start2 * self.N2], color='black')
        plt.axvline(x=self.T2[end2 * self.N2], color='black')
        f.show()





    def flattenup(self, words, p):
        curr_len = 0
        index = 0
        for i in words:
            if i == 0:
                curr_len += 1
            else:
                if p > curr_len > 0:
                    for j in range(index - curr_len, index):
                        words[j] = 1
                curr_len = 0
            index += 1

    def flattendown(self, words, q):
        curr_len = 0
        index = 0
        for i in words:
            if i == 1:
                curr_len += 1
            else:
                if q > curr_len > 0:
                    for j in range(index - curr_len, index):
                        words[j] = 0
                curr_len = 0
            index += 1



    def drawBorders(self, words, data, T, N):
        start = words.index(1)
        end = len(words) - 1 - words[::-1].index(1)

        rec = data[start * N: end * N]

        f = plt.figure()
        plt.plot(T[::], data[::], 'b')
        plt.axvline(x=T[start * N], color='r')
        plt.axvline(x=T[end * N], color='r')
        f.show()
        print(rec)
        return rec

    def furije(self,y, N):

        y_temp = np.fft.fft(y)[0:int(N / 2)] / N
        y_temp[1:] = 2 * y_temp[1:]
        FFT_y = np.abs(y_temp)
        #np.fft.fft(self.data)
        return FFT_y

    def window_function(self):
        if(self.win_fun.get() == 0):
            print("Hanning")

            self.rec = self.rec * np.hanning(len(self.rec))


        elif (self.win_fun.get() == 1):
            print("Hamming")

            self.rec = self.rec * np.hamming(len(self.rec))


        elif (self.win_fun.get() == 2):
            print("None")

            self.rec = self.rec * np.ones(len(self.rec))


    def drawNoise(self):


        f = plt.figure()
        plt.plot(self.T[::], self.data[::], 'b')
        f.show()



    def histogram(self):

        self.window_function()

        FFT_rec = self.furije(self.rec, self.N)

        self.freq = np.linspace(0,  self.N, self.N)


        self.fft_rec = np.interp(FFT_rec, (FFT_rec.min(), FFT_rec.max()), (0, 100))

        h1 = abs(self.fft_rec[0:self.N])
        h2 = np.flip(h1)
        self.magn = np.append(h1, h2)

        plt.plot(self.freq, self.magn)
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Magnitude (%)")
        plt.show()

    def spectogram(self):

        if len(self.data) == 0:
            messagebox.showinfo("Error", "File not processed")
            return
        print(self.data)
        f = plt.figure()
        plt.title("Spectrogram")
        plt.xlabel("Time (s)")
        plt.ylabel("Frequency")
        plt.specgram(self.data, Fs=self.samplerate)


        try:
            begin = self.words.index(1) - 1
            end = len(self.words) - self.words[::-1].index(1)
            plt.axvline(x=self.T[begin * self.N], color='r')
            plt.axvline(x=self.T[end * self.N], color='r')
            f.show()
        except:
            messagebox.showinfo("Error", "File is noise")




root = Tk()
root.resizable(False, False)
App(root)
root.mainloop()