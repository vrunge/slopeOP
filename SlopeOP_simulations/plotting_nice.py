import matplotlib.pyplot as plt
import os

# SET FONTS, COLORMAPS etc
# plt.rcParams["figure.figsize"] = (15,5)
plt.rcParams["image.cmap"] = "hot"
plt.rcParams["image.interpolation"] = "nearest"
plt.rcParams["font.size"] = 14

# save fig parameters
argv = {
    'dpi':600,
    'transparent':True, # (sub)plots have a transparent background
    'facecolor':'white', # The figure background
    'bbox_inches':'tight',
    'pad_inches':0}

def savefig(paht, alerts=False, **kwargs):
    """Save the current figure for printing."""
    gca = plt.gca()
    if alerts:
        assert gca.get_title()!='', "Missing figure title"
        assert gca.get_xlabel()!='', "Missing x-axis label"
        assert gca.get_ylabel()!='', "Missing y-axis label"

    plt.tight_layout()

    _, file_extension = os.path.splitext(paht)
    assert file_extension.lower() in ['.png','.svg','.pdf'], "possible extensions are .png or .svg"
    argv['format']=file_extension.lower()[1:]
    
    for k,v in kwargs.items():
        argv[k]=v
    
    plt.savefig(paht, **argv)
    