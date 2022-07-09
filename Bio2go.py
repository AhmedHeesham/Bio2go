import tkinter
import tkinter.messagebox
import customtkinter
import os
import re
import sys
from Bio.pairwise2 import format_alignment
from Bio import pairwise2
from Bio.Seq import Seq
import pylab
from Bio import Phylo
import Bio as Bio
from io import StringIO
from Bio.Phylo.NewickIO import Parser
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import matplotlib
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo

customtkinter.set_appearance_mode("System")  # Modes: "System" (standard), "Dark", "Light"
customtkinter.set_default_color_theme("blue")  # Themes: "blue" (standard), "green", "dark-blue"


class App(customtkinter.CTk):
    WIDTH = 780
    HEIGHT = 520

    def __init__(self):
        super().__init__()

        self.title("Bio2go")
        self.geometry(f"{App.WIDTH}x{App.HEIGHT}")
        self.protocol("WM_DELETE_WINDOW", self.on_closing)  # call .on_closing() when app gets closed

        # ============ create two frames ============

        # configure grid layout (2x1)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        self.frame_left = customtkinter.CTkFrame(master=self,
                                                 width=180,
                                                 corner_radius=0)
        self.frame_left.grid(row=0, column=0, sticky="nswe")

        self.frame_right = customtkinter.CTkFrame(master=self)
        self.frame_right.grid(row=0, column=1, sticky="nswe", padx=20, pady=20)

        # ============ frame_left ============

        # configure grid layout (1x11)
        self.frame_left.grid_rowconfigure(0, minsize=10)  # empty row with minsize as spacing
        self.frame_left.grid_rowconfigure(5, weight=1)  # empty row as spacing
        self.frame_left.grid_rowconfigure(8, minsize=20)  # empty row with minsize as spacing
        self.frame_left.grid_rowconfigure(11, minsize=10)  # empty row with minsize as spacing

        self.label_1 = customtkinter.CTkLabel(master=self.frame_left,
                                              text="Tools",
                                              text_font=("Roboto Medium", -16))  # font name and size in px
        self.label_1.grid(row=0, column=0, pady=10, padx=10)

        self.button_1 = customtkinter.CTkButton(master=self.frame_left,
                                                text="fastq to fasta",
                                                command=self.FastqToFasta)
        self.button_1.grid(row=1, column=0, pady=10, padx=20)

        self.button_2 = customtkinter.CTkButton(master=self.frame_left,
                                                text="Multiple sequence ",
                                                command=self.MulitpleSequanceAlignment)
        self.button_2.grid(row=2, column=0, pady=10, padx=20)

        self.button_3 = customtkinter.CTkButton(master=self.frame_left,
                                                text="Alignment",
                                                command=self.Alignment)
        self.button_3.grid(row=3, column=0, pady=10, padx=20)
        self.button_4 = customtkinter.CTkButton(master=self.frame_left,
                                                text="Validation",
                                                command=self.rawan)
        self.button_4.grid(row=4, column=0, pady=10, padx=20)
        self.button_5 = customtkinter.CTkButton(master=self.frame_left,
                                                text="Creating a tree",
                                                command=self.Tree)
        self.button_5.grid(row=5, column=0, pady=10, padx=20)

        self.button_6 = customtkinter.CTkButton(master=self.frame_left,
                                                text="Some statial operations",
                                                command=self.statialOperations)
        self.button_6.grid(row=6, column=0, pady=10, padx=20)
        self.button_8 = customtkinter.CTkButton(master=self.frame_left,
                                                text="Blast",
                                                command=self.zaharaa)
        self.button_8.grid(row=7, column=0, pady=10, padx=20)
        self.button_7 = customtkinter.CTkButton(master=self.frame_right,
                                                text="Automata",
                                                border_width=2,  # <- custom border_width
                                                fg_color=None,  # <- no fg_color
                                                command=self.allOnOne)
        self.button_7.grid(row=8, column=2, columnspan=1, pady=20, padx=20, sticky="we")

        self.label_mode = customtkinter.CTkLabel(master=self.frame_left, text="Appearance Mode:")
        self.label_mode.grid(row=9, column=0, pady=0, padx=20, sticky="w")

        self.optionmenu_1 = customtkinter.CTkOptionMenu(master=self.frame_left,
                                                        values=["Light", "Dark", "System"],
                                                        command=self.change_appearance_mode)
        self.optionmenu_1.grid(row=10, column=0, pady=10, padx=20, sticky="w")

        # ============ frame_right ============

        # configure grid layout (3x7)
        self.frame_right.rowconfigure((0, 1, 2, 3), weight=1)
        self.frame_right.rowconfigure(7, weight=10)
        self.frame_right.columnconfigure((0, 1), weight=1)
        self.frame_right.columnconfigure(2, weight=0)

        self.frame_info = customtkinter.CTkFrame(master=self.frame_right)
        self.frame_info.grid(row=0, column=0, columnspan=2, rowspan=4, pady=20, padx=20, sticky="nsew")

        # ============ frame_info ============

        # configure grid layout (1x1)
        self.frame_info.rowconfigure(0, weight=1)
        self.frame_info.columnconfigure(0, weight=1)

        self.label_info_1 = customtkinter.CTkLabel(master=self.frame_info,
                                                   text="it is Bio2go\n" +
                                                        " \n",
                                                   height=100,
                                                   fg_color=("white", "gray38"),  # <- custom tuple-color
                                                   justify=tkinter.LEFT)
        self.label_info_1.grid(column=0, row=0, sticky="nwe", padx=15, pady=15)

        # ============ frame_right ============

        self.radio_var = tkinter.IntVar(value=0)

        self.label_radio_group = customtkinter.CTkLabel(master=self.frame_right,
                                                        text="Future work:")
        self.label_radio_group.grid(row=0, column=2, columnspan=1, pady=20, padx=10, sticky="")

        # set default values
        self.optionmenu_1.set("Dark")

    def button_event(self):
        print("Button pressed")

    def MulitpleSequanceAlignment(self):
        dialogOptions = customtkinter.CTkInputDialog(master=None, text="1-retreve data from ncbi\n2-you have a file:", title="MSA")
        options = dialogOptions.get_input()
        if options == "1":
            lst = []
            dialog = customtkinter.CTkInputDialog(master=None, text="Enter number of elements:", title="Number of elements")
            n = int(dialog.get_input())
            n=int(n)
            for i in range(0, n):
                dialog1 = customtkinter.CTkInputDialog(master=None, text="Enter the accession number:", title="Accession number")
                ele = dialog1.get_input()
                lst.append(ele)
            y = ','.join(lst)
            print(y)
            InputFileName = customtkinter.CTkInputDialog(master=None, text="Enter the retrieved file name:", title="Accession number")
            nameOfFile = InputFileName.get_input()
            command = "efetch -db nucleotide -format fasta -id {} >{}.fasta".format(y,nameOfFile)
            os.system(command)

        elif options =="2":
            dialogFile = customtkinter.CTkInputDialog(master=None, text="Enter file name:", title="File name")
            nameOfFile = dialogFile.get_input()
        dialog2 = customtkinter.CTkInputDialog(master=None, text="Enter format\'for example clw\':", title="Format")
        e = dialog2.get_input()
        nameOfOutFile = customtkinter.CTkInputDialog(master=None, text="Enter the out file name:", title="Format")
        nameOfOutFile = nameOfOutFile.get_input()
        command2 = "muscle -in {}.fasta -out {}.{} -{}".format(nameOfFile,nameOfOutFile,e,e)
        os.system(command2)
        window = customtkinter.CTkToplevel(self)
        window.geometry("400x200")
        window.title("output")
        label = customtkinter.CTkLabel(master=window, text="Your file created\nclose this window <3")
        label.place(relx=0.5, rely=0.5, anchor=tkinter.CENTER)


    def FastqToFasta(self):
        dialog = customtkinter.CTkInputDialog(master=None, text="Pleaser Enter your file name:", title="FastqToFasta")
        inFile = dialog.get_input()
        dialog1 = customtkinter.CTkInputDialog(master=None, text="Pleaser Enter your out file name:", title="FastqToFasta")
        outFile = dialog1.get_input()
        command1 = 'sed -n \'1~4s/^@/>/p;2~4p\' {}.fastq > {}.fasta'.format(inFile,outFile)
        os.system(command1)
        window = customtkinter.CTkToplevel(self)
        window.geometry("400x200")
        window.title("output")
        label = customtkinter.CTkLabel(master=window, text="Your file created\nclose this window <3")
        label.place(relx=0.5, rely=0.5, anchor=tkinter.CENTER)

    '''
    def FastqToFasta(self):
        x = input("Please Enter your file name: ")
        command1 = 'sed -n \'1~4s/^@/>/p;2~4p\' {}.fastq > OUTFILE.fasta'.format(
            x)
        os.system(command1)
        print("your file created <3")
    '''

    def Alignment(self):
        dialogOp = customtkinter.CTkInputDialog(master=None, text="1-golbal alignment\n 2-local alignment:", title="options")
        op = dialogOp.get_input()
        if op == "1":
            dialog1 = customtkinter.CTkInputDialog(master=None, text="Enter seq1:", title="seq1")
            seq1 = Seq(dialog1.get_input())
            dialog2 = customtkinter.CTkInputDialog(master=None, text="Enter seq2:", title="seq2")
            seq2 = Seq(dialog2.get_input())
            alignments = pairwise2.align.globalxx(seq1, seq2)
            for alignment in alignments:
                print(format_alignment(*alignment))
        if op == "2":
            dialog1 = customtkinter.CTkInputDialog(master=None, text="Enter seq1:", title="seq1")
            seq1 = Seq(dialog1.get_input())
            dialog2 = customtkinter.CTkInputDialog(master=None, text="Enter seq2:", title="seq2")
            seq2 = Seq(dialog2.get_input())
            alignments = pairwise2.align.localxx(seq1, seq2)
            for alignment in alignments:
                print(format_alignment(*alignment))        
    '''
    def align_global(self):
        x = input("Enter seq1: ")
        y = input("Enter seq2: ")
        seq1 = Seq(x)
        seq2 = Seq(y)

        alignments = pairwise2.align.globalxx(seq1, seq2)
        for alignment in alignments:
            print(format_alignment(*alignment))
    '''
    def allOnOne(self):
        dialogOptions = customtkinter.CTkInputDialog(master=None, text="1-retreve data from ncbi\n2-you have a file:", title="MSA")
        options = dialogOptions.get_input()
        if options == "1":
            lst = []
            dialog = customtkinter.CTkInputDialog(master=None, text="Enter number of elements:", title="Number of elements")
            n = int(dialog.get_input())
            for i in range(0, n):
                dialog1 = customtkinter.CTkInputDialog(master=None, text="Enter the accession number:", title="Accession number")
                ele = dialog1.get_input()
                lst.append(ele)
            y = ','.join(lst)
            print(y)
            InputFileName = customtkinter.CTkInputDialog(master=None, text="Enter the retrieved file name:", title="Accession number")
            nameOfFile = InputFileName.get_input()
            command = "efetch -db nucleotide -format fasta -id {} >{}.fasta".format(y,nameOfFile)
            os.system(command)

        elif options =="2":
            dialogFile = customtkinter.CTkInputDialog(master=None, text="Enter file name:", title="File name")
            nameOfFile = dialogFile.get_input()
        nameOfOutFile = customtkinter.CTkInputDialog(master=None, text="Enter the out file name:", title="Format")
        nameOfOutFile = nameOfOutFile.get_input()
        command2 = "muscle -in {}.fasta -out {}.clw -clw".format(nameOfFile,nameOfOutFile)
        os.system(command2)
        with open("{}.clw".format(nameOfOutFile), "r") as clw:
            alignment = AlignIO.read(clw, "clustal")
        # print(type(alignment))
        calculator = DistanceCalculator('identity')
        distance_matrix = calculator.get_distance(alignment)
        # print(distance_matrix)
        constructor = DistanceTreeConstructor(calculator)
        turtle_tree = constructor.build_tree(alignment)
        turtle_tree.rooted = True
        # print(turtle_tree)
        Phylo.write(turtle_tree, "nadeen.xml", "phyloxml")
        # fig = Phylo.draw(turtle_tree)
        fig = plt.figure(figsize=(13, 5), dpi=100)  # create figure & set the size
        matplotlib.rc('font', size=12)  # fontsize of the leaf and node labels
        matplotlib.rc('xtick', labelsize=10)  # fontsize of the tick labels
        matplotlib.rc('ytick', labelsize=10)  # fontsize of the tick labels
        # turtle_tree.ladderize()
        axes = fig.add_subplot(1, 1, 1)
        Phylo.draw(turtle_tree, axes=axes)

        fig.savefig("nadeen")
        

    def aotomated(self):
        window = customtkinter.CTkToplevel(self)
        window.geometry("586x390")
        window.title("Aotomata")
        window.grid_columnconfigure(1, weight=1)
        window.grid_rowconfigure(0, weight=1)
        window.frame_left = customtkinter.CTkFrame(master=window,
                                                   width=293,
                                                   corner_radius=0)
        window.frame_left.grid(row=0, column=0, sticky="nswe")

        window.frame_right = customtkinter.CTkFrame(master=window)
        window.frame_right.grid(row=0, column=1, sticky="nswe", padx=2, pady=0)
    
        '''
        lst = []
        dialog = customtkinter.CTkInputDialog(master=None, text="Enter number of elements:", title="Number of elements")
        n = int(dialog.get_input())
        for i in range(0, n):
            dialog1 = customtkinter.CTkInputDialog(master=None, text="Enter the accession number:", title="Accession number")
            ele = dialog1.get_input()
            lst.append(ele)
        y = ','.join(lst)
        print(y)
        InputFileName = customtkinter.CTkInputDialog(master=None, text="Enter the retrieved file name:", title="Accession number")
        nameOfFile = InputFileName.get_input()
        command = "efetch -db nucleotide -format fasta -id {} >{}.fasta".format(y,nameOfFile)
        dialog2 = customtkinter.CTkInputDialog(master=None, text="Enter format\'for example clw\':", title="Format")
        e = dialog2.get_input()
        nameOfOutFile = customtkinter.CTkInputDialog(master=None, text="Enter the out file name:", title="Format")
        nameOfOutFile = nameOfOutFile.get_input()
        command2 = "muscle -in {}.fasta -out {}.{} -{}".format(nameOfFile,nameOfOutFile,e,e)
        os.system(command)
        os.system(command2)
        window = customtkinter.CTkToplevel(self)
        window.geometry("400x200")

        label = customtkinter.CTkLabel(master=window, text="Your file created\nclose this window <3")
        label.place(relx=0.5, rely=0.5, anchor=tkinter.CENTER)
'''

    def statialOperations(self):
        op = customtkinter.CTkInputDialog(master=None,
                                          text="what is the statistical operation do you want to do\n1- calculate nucleotide\n2-calculate intrions and exons\n Enter the number to choose:",
                                          title="operations")
        op = op.get_input()
        if op == "1":
            x = customtkinter.CTkInputDialog(master=None, text="enter the nucleotide you want to calculate ACGT:",
                                             title="nucleotide")
            x = x.get_input()
            y = customtkinter.CTkInputDialog(master=None, text="enter the file name:", title="file name")
            y = y.get_input()
            x.upper()
            if x == "A" or "C" or "G" or "T":
                command = "grep -o -i {} {} | wc -l".format(x, y)
                os.system(command)
            else:
                print("enter a valid nucleotide")

        elif op == "2":
            y = input("enter the file name:")
            command = "gunzip -c {}| cut -f 3 | sort | uniq -c".format(y)
            os.system(command)

    def Tree(self):
        dialog = customtkinter.CTkInputDialog(master=None, text="Pleaser Enter your file name :", title="Creat a tree")
        with open("{}.clw".format(dialog.get_input()), "r") as clw:
            alignment = AlignIO.read(clw, "clustal")
        # print(type(alignment))
        calculator = DistanceCalculator('identity')
        distance_matrix = calculator.get_distance(alignment)
        # print(distance_matrix)
        constructor = DistanceTreeConstructor(calculator)
        turtle_tree = constructor.build_tree(alignment)
        turtle_tree.rooted = True
        # print(turtle_tree)
        Phylo.write(turtle_tree, "nadeen.xml", "phyloxml")
        # fig = Phylo.draw(turtle_tree)
        fig = plt.figure(figsize=(13, 5), dpi=100)  # create figure & set the size
        matplotlib.rc('font', size=12)  # fontsize of the leaf and node labels
        matplotlib.rc('xtick', labelsize=10)  # fontsize of the tick labels
        matplotlib.rc('ytick', labelsize=10)  # fontsize of the tick labels
        # turtle_tree.ladderize()
        axes = fig.add_subplot(1, 1, 1)
        Phylo.draw(turtle_tree, axes=axes)

        fig.savefig("nadeen")
        '''
        dialog = customtkinter.CTkInputDialog(master=None, text="Pleaser Enter your file name :", title="Creat a tree")
        with open("{}.clw".format(dialog.get_input()), "r") as clw:
            alignment = AlignIO.read(clw, "clustal")
        # print(type(alignment))
        calculator = DistanceCalculator('identity')
        distance_matrix = calculator.get_distance(alignment)
        # print(distance_matrix)
        constructor = DistanceTreeConstructor(calculator)
        turtle_tree = constructor.build_tree(alignment)
        turtle_tree.rooted = True
        # print(turtle_tree)
        Phylo.write(turtle_tree, "turtle_tree.xml", "phyloxml")
        # fig = Phylo.draw(turtle_tree)
        fig = plt.figure(figsize=(13, 5), dpi=100)  # create figure & set the size
        matplotlib.rc('font', size=12)  # fontsize of the leaf and node labels
        matplotlib.rc('xtick', labelsize=10)  # fontsize of the tick labels
        matplotlib.rc('ytick', labelsize=10)  # fontsize of the tick labels
        # turtle_tree.ladderize()
        axes = fig.add_subplot(1, 1, 1)
        Phylo.draw(turtle_tree, axes=axes)

        fig.savefig("turtles_cladogram")
'''
    def change_appearance_mode(self, new_appearance_mode):
        customtkinter.set_appearance_mode(new_appearance_mode)

    def on_closing(self, event=0):
        self.destroy()
    def zaharaa(self):
            dialogFileName = customtkinter.CTkInputDialog(master=None, text="Enter the file name of database:", title="File name")
            x = dialoFileName.get_input()
            command1 = "makeblastdb -in {} -dbtype nucl".format(x)
            lst = []
            dialogelemet = customtkinter.CTkInputDialog(master=None, text="Enter number of elements:", title="Number of elements")
            n = dialogelemet.get_input()
            for i in range(0, n):
                ele = input()
                lst.append(ele)
            y=','.join(lst)
            print(y)
            for j in lst:
                command2 = "blastn -query {}.fa -db {} -out {}.txt".format(j, x, j)
            os.system(command1)
            os.system(command2) 
            window = customtkinter.CTkToplevel(self)
            window.geometry("400x200")
            window.title("output")
            label = customtkinter.CTkLabel(master=window, text="html file is created\nclose this window <3")
            label.place(relx=0.5, rely=0.5, anchor=tkinter.CENTER)
    def rawan(self):
        dialogOptions = customtkinter.CTkInputDialog(master=None, text="1- check the quality\n2-validat fastq file :", title="options")
        op = dialogOptions.get_input() 
        if op=='1':
            dialoFileName = customtkinter.CTkInputDialog(master=None, text="Enter file name:", title="Quality")
            fileName = dialoFileName.get_input()
            command1 = 'zcat {}.gz | fastqc stdin'.format(fileName)
            os.system(command1)
            window = customtkinter.CTkToplevel(self)
            window.geometry("400x200")
            window.title("output")
            label = customtkinter.CTkLabel(master=window, text="html file is created\nclose this window <3")
            label.place(relx=0.5, rely=0.5, anchor=tkinter.CENTER)
        elif op=='2':    
            class SeqRead:
                # Initiate object in this class with its attributes
                def __init__(self, read_id, seq, plus, qualit):
                    self.id = read_id
                    self.seq = seq
                    self.plus = plus
                    self.quality = qualit

                # start analysis with checking header line starts with "@"
                def header_control(self, report, line_no):
                    if not self.id.startswith('@'):
                        report.write('invalid header.{}\n'.format(line_no - 3))

                # checking sequence line if it contains any other character other than  "ATCGNatcgn"
                def seq_errors(self, report, line_no):
                    if not re.match('[ATCGNatcgn]+$', self.seq):
                        report.write('Invalid seq.{}\n'.format(line_no - 2))

                # checking plus line if it contains + or not
                def no_plus(self, report, line_no):
                    if self.plus != '+\n':
                        report.write('Invalid plus_line.{}\n'.format(line_no - 1))

                # checking quality line if its length of characters  = length of sequence line 
                def qualine_errors(self, report, line_no):
                    if len(self.quality) != len(self.seq):
                        report.write('Invalid quality.{}\n'.format(line_no))

            # open &read compressed file
            def read_gz(file_path):
                if file_path.endswith('.gz'):
                    file = gzip.open(file_path)
                else:
                    file = open(file_path, 'r')
                return file

            # specify the path of your file

            # input file 'fastq'
            #file_path = "FASTQ_example.fastq"
            dialog = customtkinter.CTkInputDialog(master=None, text="Enter your file Name :", title="File name")
            file_path = dialog.get_input()
            # min_len= int(int(sys.argv[0]))

            # Open the file with keyword 'with', with open let file to be closed after iteration
            f = read_gz(file_path)
            # report_file= sys.argv[2]
            r_name = file_path + "validate_fastaq"

            counter = 0
            with open(r_name, 'w') as report:
                # A loop that reads each 4 lines in the file as read object and puts them in a list of read objects
                while True:
                    lines = []
                    for i in range(4):
                        line = f.readline()
                        lines.append(line)
                        counter += 1
                    if not line:
                        break
                    # Create a read object named SeqRead
                    read = SeqRead(read_id=lines[0], seq=lines[1], plus=lines[2], qualit=lines[3])
                    read.header_control(report, counter)
                    read.seq_errors(report, counter)
                    read.no_plus(report, counter)
                    read.qualine_errors(report, counter)
            window = customtkinter.CTkToplevel(self)
            window.geometry("400x200")
            window.title("output")
            label = customtkinter.CTkLabel(master=window, text="A file created \'if it is embty your file is a fastq file\'")
            label.place(relx=0.5, rely=0.5, anchor=tkinter.CENTER)

if __name__ == "__main__":
    app = App()
    app.mainloop()
