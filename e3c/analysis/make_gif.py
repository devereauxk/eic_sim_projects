import os, sys
import fitz
from PIL import Image
import glob

# Get the PDF directory from command-line argument
if len(sys.argv) < 2:
    print("Usage: python script_name.py pdf_directory")
    sys.exit(1)

def animate_dir(pdf_directory):
    print("processing " + pdf_directory + " ...")

    # Temporary directory to store PNG files
    png_temp_directory = "./png_temp"

    # List of X values
    x_values = range(40, 51)

    # Create the temporary directory if it doesn't exist
    if not os.path.exists(png_temp_directory):
        os.makedirs(png_temp_directory)

    # Convert PDF files to PNG files
    frames = []
    for X in x_values:
        print(X)
        pdf_filename = f"h2d_jet_e3c_xi_phi_{X}_surface.pdf"
        pdf_path = os.path.join(pdf_directory, pdf_filename)
        
        png_filename = f"{X}.png"
        png_path = os.path.join(png_temp_directory, png_filename)
        
        doc = fitz.open(pdf_path)  # open document
        for i, page in enumerate(doc):
            pix = page.get_pixmap(dpi=300)  # render page to an image
            pix.save(png_path)
            break

        new_frame = Image.open(png_path)
        frames.append(new_frame)


    # Save into a GIF file that loops forever
    frames[0].save(os.path.join(pdf_directory,'surface.gif'), format='GIF',
                    append_images=frames[1:],
                    save_all=True,
                    duration=200, loop=0)

    # Clean up: Remove temporary PNG files
    for X in x_values:
        png_filename = os.path.join(png_temp_directory, f"{X}.png")
        os.remove(png_filename)

    # Remove temporary directory
    os.rmdir(png_temp_directory)
    

all_dirs = [
    "./ep_10_100_K0_pow025/", "./ep_10_100_K0_pow05/", "./ep_10_100_K0_pow1/", "./ep_10_100_K0_pow15/",
    "./eAu_10_100_K4_density_pow025/", "./eAu_10_100_K4_density_pow05/", "./eAu_10_100_K4_density_pow1/", "./eAu_10_100_K4_density_pow15/", 
    "./eU_10_100_K4_density_pow025/", "./eU_10_100_K4_density_pow05/", "./eU_10_100_K4_density_pow1/", "./eU_10_100_K4_density_pow15/"
]

pdf_directory = sys.argv[1]

if pdf_directory == "all":
    for dir in all_dirs:
        animate_dir(dir)
else:
    animate_dir(pdf_directory)