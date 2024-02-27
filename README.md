# Brainviz

A python script to generate text files with XYZ coordinates from various data formats.   

## Usage

- Install Python and run something like `py export_brain_viz.py` or `python export_brain_viz.py` depending on your main Python CLI call.  
- Edit file according to your needs. The *data* directory contains sources, while *output* contains your exported text files.  
- The *data* directory is left empty, for you to fill it with relevant brain data files.  
- Adapt the script by commenting or uncommenting revelant code sections.
- The **dn** parameters is used to generate sparser outputs (a **dn** of 2 will generate a file that is 2³ == 8 times lighter), and the **threshold** parameter filters points with a minimum value

Supported formats are currently **.raw**, **.nii**, **.trk** and **.tif**. The code was set up to visualize brain masks and fiber segments during a dome show by Fabien Chauveau and Joshua Gobé, held on *March 21st 2024*.   
 
