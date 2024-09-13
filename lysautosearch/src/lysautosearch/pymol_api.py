#TKINKER IMPORTS
import functools
import threading
import tkinter
from collections import defaultdict
from tkinter import *
import numpy as np
import pandas as pd
import pymol
from pymol.cgo import *
from pymol.vfont import plain
import sys
from typing import Union
#Readline
try:
  import readline #Linux
except ImportError:
  import pyreadline as readline #Windows
def call_pymol(sname:str,
               pdb_file:str,
               mapped_positions_file:str,
               results_dir:str,
               vars_dict
               ):
    '''Visualization program

    :params
    lab_x: Positions of the label in the x-axis
    lab_y: Positions of the label in the y-axis
    lab_z: Positions of the label in the z-axis

    '''

    if vars_dict is not None:
        pass
    else:
        vars_dict = {"angle_x":-50, "angle_y":210, "angle_z":20, "label_x":70, "label_y":50, "label_z":80, "zoom":-5, "shape":"cartoon"}

    angle_x = vars_dict["angle_x"]
    angle_y = vars_dict["angle_y"]
    angle_z = vars_dict["angle_z"]
    lab_x = vars_dict["label_x"]
    lab_y = vars_dict["label_y"]
    lab_z = vars_dict["label_z"]
    zoom = vars_dict["zoom"]
    shape = vars_dict["shape"]

    #LAUNCH PYMOL
    #readline.parse_and_bind('tab: complete')  # Set up path to pymol environment (in case is not installed)
    # pymol launching: quiet (-q), without GUI (-c) and with arguments from command line
    #pymol.pymol_argv = ['pymol', '-q'] + sys.argv[1:]
    pymol.finish_launching()


    # Load Structures
    pymol.cmd.load(pdb_file, sname)
    pymol.cmd.disable("all")  # toggles off the display of all currently visible representations of an object. It is the equivalent of deselecting the object
    pymol.cmd.enable(sname)


    def colour_by_selection(selection="all",
                            Selected="orange",
                            Not='grey50',
                            Domain='lime',
                            Selected_and_Domain='magenta',
                            ):

        colors = {
            'Selected': Selected,
            'Not': Not,
            'Domain': Domain,
            'Selected_and_Domain': Selected_and_Domain
        }

        #BACKGROUND & SHAPES
        cmd.bg_color('white')
        cmd.show_as(shape, 'all')
        cmd.color('gray', 'all')
        #ROTATION
        cmd.rotate([1, 0, 0], angle=angle_x, selection="all")
        cmd.rotate([0, 1, 0], angle = angle_y,selection = "all") # Commands to rotate the structures to visualize some specific side of the protein [x,y,z]
        cmd.rotate([1, 0, 0], angle=angle_z, selection="all")
        #ELIMINATING CHAINS
        # Eliminating chains in the structure if desired
        # cmd.select('chainA', 'chain A')
        # cmd.remove('chain A')

        #LEGEND
        ###The text that appears in the image, change placement accordingly
        cgo = []
        axes = [[5.0, 0.0, 0.0], [0.0, 5.0, 0.0],[0.0, 0.0, 5.0]]  # Change the values if the protein does not quite fall into place


        cyl_text(cgo, plain, [float(lab_x),float(lab_y), float(lab_z)], '%s' % (sname.split('.')[0]), radius=0.6, color=[0.0, 0.0, 0.0],axes=axes)  # x=60 for RIOK2, x=40 and z=60 for ROS1
        cyl_text(cgo, plain, [float(lab_x), float(lab_y)-10.0, float(lab_z)], 'Positively Selected', radius=0.6, color=[1.0, 0.5, 0.0], axes=axes)
        cyl_text(cgo, plain, [float(lab_x), float(lab_y)-20.0, float(lab_z)], 'Not selected', radius=0.6, color=[0.5, 0.5, 0.5], axes=axes)
        cyl_text(cgo, plain, [float(lab_x), float(lab_y)-30.0, float(lab_z)], 'Functional Domain', radius=0.6, color=[0.5, 1.0, 0.5], axes=axes)
        cyl_text(cgo, plain, [float(lab_x), float(lab_y)-40.0, float(lab_z)], 'Both', radius=0.6, color=[1.0, 0.0, 1.0], axes=axes)

        cmd.set("cgo_line_radius", 0.03)  # 0.03
        cmd.load_cgo(cgo, 'txt')
        #ZOOM
        cmd.zoom("all", zoom)  # Higher and positive values zoom  out, it accepts negative values

        #READ PREVIOUSLY CREATED DATAFRAME:

        mapped_positions = pd.read_csv(mapped_positions_file,sep="\t")
        #Option A: best alternative
        mapped_positions['PDB_Position'] = mapped_positions['PDB_Position'].astype(np.float64) #Need to convert to float to use isfinite
        mapped_positions = mapped_positions[np.isfinite(mapped_positions['PDB_Position'])]  # Remove the Residues that got 'nan' value in their equivalent positions
        position_phenotype_dict = pd.Series(mapped_positions.Label.values, index=mapped_positions.PDB_Position).to_dict()
        # Colour the residues in the protein according to their label in the dataframe
        for key, value in position_phenotype_dict.items():
            #print(int(key), value, colors[value])
            cmd.color(colors[value], 'resi %s' % int(key))  # --------->If it does not work (no colour shown) probably it has to do with the Residues ID being wrong
            ###LABEL
            if value == 'Selected_and_Domain': #Label the alpha carbon from positions selected and in the domain
                #print(key)
                #cmd.select('Both','resi %s' % int(key)) #Create a selection
                #cmd.label('Both and n. CA', '" %s %s" % (resn,resi)')
                cmd.label('resi %s and n. CA' % int(key), '" %s %s" % (resn,resi)')
                #cmd.label('resi %s' % int(key), '" %s %s" % (resn,resi)')


    cmd.extend("Colour_by_Selection", colour_by_selection)
    print("Structure will be at %s" % (results_dir))
    colour_by_selection(sname)

    pymol.cmd.png(f"{results_dir}/{sname.split('.')[0]}_{shape}",ray=1)
    #pymol.cmd.png(results_dir + "/%s_{}" % (sname.split('.')[0],shape),ray=1)

def transform(variable):
    if not isinstance(variable,str):
        variable = eval(variable)
    return variable


class MyFrame(Frame):
    def __init__(self,den,vars_dict,dict_keys):
        Frame.__init__(self,den)
        self.vars_dict = vars_dict
        self.dict_keys = dict_keys
        self.master.rowconfigure(10, weight=1)
        self.master.columnconfigure(4, weight=0) #weight 0 is supposed to avoid column spread
        #self.master.maxsize(300, 300)

        #Rotate protein (x,y,z entries)
        self.angle_x_label = Label(den,text='Rotation axis x')
        self.angle_y_label = Label(den,text='Rotation axis y')
        self.angle_z_label = Label(den,text='Rotation axis z')
        self.angle_x = Entry(den,width=5)
        self.angle_y = Entry(den,width=5)
        self.angle_z = Entry(den,width=5)

        #Label placement (x,y,z entries)
        self.label_label = Label(den, text='Label Placement')
        self.label_x = Entry(den,width=5)
        self.label_y = Entry(den,width=5)
        self.label_z = Entry(den,width=5)
        #Zoom
        self.zoom_label = Label(den,text='Zoom')
        self.zoom_entry= Entry(den,width=5)
        #Residues shape
        self.shape_label = Label(den, text='Shape')
        self.shape_entry = Entry(den, width=8)
        #Default values shown to user

        self.angle_x.insert(0,20)
        self.angle_y.insert(0,0)
        self.angle_z.insert(0,0)
        self.label_x.insert(0,70)
        self.label_y.insert(0,50)
        self.label_z.insert(0,80)
        self.zoom_entry.insert(0,10)
        self.shape_entry.insert(0,'cartoon')

        ##RUN
        self.run = Button(den, text="RUN", command=lambda: [f() for f in [self.button_action, self.quit,self.destroy]],width=40, background='green')

        # self.master.columnconfigure(0, weight=1)  # First column can expand
        # self.master.columnconfigure(1, weight=0)  # Second column won't expand
        # self.master.columnconfigure(2, weight=0)  # Third column can expand
        # self.master.columnconfigure(3, weight=0)  # Third column can expand

        self.angle_x_label.grid(row=0, column=0)
        self.angle_y_label.grid(row=1, column=0)
        self.angle_z_label.grid(row=2, column=0)
        self.angle_x.grid(row=0, column=1)
        self.angle_y.grid(row=1, column=1)
        self.angle_z.grid(row=2, column=1)
        self.label_label.grid(row=3, column=0)
        self.label_x.grid(row=3, column=1,sticky="ew")
        self.label_y.grid(row=3, column=2,sticky="ew")
        self.label_z.grid(row=3, column=3,sticky="ew")
        self.zoom_label.grid(row=4, column=0)
        self.zoom_entry.grid(row=4, column=1)
        self.shape_label.grid(row=5, column=0)
        self.shape_entry.grid(row=5, column=1)
        self.run.grid(row=8, column=1)

    def button_action(self): #Insert all the default values here
        pymol.cmd.reinitialize('everything')

        for key,val in self.dict_keys.items():
            if hasattr(self,val[0]):
                getval = getattr(self,val[0]).get()
                setattr(self, val[1], getval)
            else:
                setattr(self, val[1], val[2])
        #self.start_pymol_thread()

    # def start_pymol_thread(self): #TODO: Make work to have tkinter and pymol Gui working together
    #     pymol_thread = threading.Thread(target=call_pymol,args=(pdb_name,
    #                pdb_file,
    #                args.lysautosearch_pymol_dataframe,
    #                args.results_dir,
    #                vars_dict))
    #     pymol_thread.daemon = True  # Ensure the thread exits when the main program exits
    #     pymol_thread.start()


def call_pymol_gui(pdb_name,pdb_file,args):

    dict_keys = {"angle_x":["angle_x","Angle_x",0],
                 "angle_y":["angle_y","Angle_y",0],
                 "angle_z":["angle_z","Angle_z",0],
                 "label_x":["label_x","lab_x",70],
                 "label_y":["label_y","lab_y",50],
                 "label_z":["label_z","lab_z",80],
                 "zoom":["zoom_entry","zoom",10],
                 "shape":["shape_entry","shape","cartoon"]}

    tmp_dict = defaultdict(list, {k: [] for k in dict_keys.keys()})

    for key,val in tmp_dict.items():
        #if val : #if list is not empty
        vals = dict_keys[key]
        tmp_dict[key] = [vals[2]] #extract the inserted value

    while True:
        # # pymol.cmd.quit()
        den = tkinter.Tk()
        den.title("Link Your Sites GUI")
        #den.geometry("500x200")

        prompt = MyFrame(den,tmp_dict,dict_keys)


        for key,val in tmp_dict.items():
            name = dict_keys[key][0]
            if val:
                getattr(prompt, name).delete(0, "end") #optional
                getattr(prompt,name).insert(0, val)

        den.mainloop()
        #Highlight: Store the input variables -> Temporary memory
        for key,val in tmp_dict.items():
            name = dict_keys[key][1]
            if val:
                del val[:]
            getval = transform(getattr(prompt,name))
            val.append(getval)

        den.destroy()  # now is working (destroying the root)
        vars_dict = {key:val[0] for key,val in tmp_dict.items()}

        call_pymol(pdb_name,
                   pdb_file,
                   args.lysautosearch_pymol_dataframe,
                   args.results_dir,
                   vars_dict)


