import os, sys
import fitz
from PIL import Image
import glob

# Get the PDF directory from command-line argument
if len(sys.argv) < 4:
    print("Usage: python script_name.py directory gif/grid fixed/not fixed")
    sys.exit(1)

def animate_dir(pdf_directory, fixed_option="fixed"):
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
    frames[0].save('./'+pdf_directory[2:len(pdf_directory)-1]+'surface.gif', format='GIF',
                    append_images=frames[1:],
                    save_all=True,
                    duration=300, loop=0)

    # Clean up: Remove temporary PNG files
    for X in x_values:
        png_filename = os.path.join(png_temp_directory, f"{X}.png")
        os.remove(png_filename)

    # Remove temporary directory
    os.rmdir(png_temp_directory)

def animate_grid(all_dirs, fixed_option="fixed"):
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
        print("processing " + str(X) + " ...")

        images = []
        diri = 0
        for pdf_directory in all_dirs:

            pdf_path = os.path.join(pdf_directory, f"h2d_jet_e3c_xi_phi_{X}_surface.pdf")
            if fixed_option == "fixed": pdf_path = os.path.join(pdf_directory, f"h2d_jet_e3c_xi_phi_{X}_surface_zfixed.pdf")
            
            png_path = os.path.join(png_temp_directory, f"{diri}.png")

            doc = fitz.open(pdf_path)  # open document
            for i, page in enumerate(doc):
                pix = page.get_pixmap(dpi=100)  # render page to an image
                pix.save(png_path)
                break
            
            img = Image.open(png_path)
            images.append(img)
            diri += 1

        # Determine the size of individual images
        image_width, image_height = images[0].size

        # Calculate the size of the final grid image
        grid_width = image_width * 4
        grid_height = image_height * 3

        # Create a new blank image for the grid
        grid_image = Image.new('RGB', (grid_width, grid_height))

        # Paste images onto the grid
        for row in range(3):
            for col in range(4):
                index = row * 4 + col
                if index < len(images):
                    grid_image.paste(images[index], (col * image_width, row * image_height))

        # Save the grid image as a PNG file
        png_path = os.path.join(png_temp_directory, f'image_grid_{X}.png')
        grid_image.save(png_path)

        new_frame = Image.open(png_path)
        frames.append(new_frame)

    # Save into a GIF file that loops forever
    output_name = './surface.gif'
    if fixed_option == "fixed": output_name = './surface_zfixed.gif'
    frames[0].save(output_name, format='GIF',
                    append_images=frames[1:],
                    save_all=True,
                    duration=700, loop=0)

    # Clean up: Remove temporary PNG files
    for idiri in range(12):
        png_filename = os.path.join(png_temp_directory, f"{idiri}.png")
        os.remove(png_filename)

    # Remove temporary directory
    os.rmdir(png_temp_directory)
    
def RL_row(pdf_directory, fixed_option="fixed"):
    print("processing " + pdf_directory + " ...")

    # Temporary directory to store PNG files
    png_temp_directory = "./png_temp"

    # List of X values
    x_values = [41, 43, 45, 47, 49, 50]

    # Create the temporary directory if it doesn't exist
    if not os.path.exists(png_temp_directory):
        os.makedirs(png_temp_directory)

    # Convert PDF files to PNG files
    images = []
    diri = 0
    for X in x_values:
        print(X)
        
        pdf_path = os.path.join(pdf_directory, f"h2d_jet_e3c_xi_phi_{X}_surface.pdf")
        if fixed_option == "fixed": pdf_path = os.path.join(pdf_directory, f"h2d_jet_e3c_xi_phi_{X}_surface_zfixed.pdf")

        png_path = os.path.join(png_temp_directory, f"{diri}.png")
        
        doc = fitz.open(pdf_path)  # open document
        for i, page in enumerate(doc):
            pix = page.get_pixmap(dpi=300)  # render page to an image
            pix.save(png_path)
            break

        new_frame = Image.open(png_path)
        images.append(new_frame)
        diri += 1

    # Determine the size of individual images
    image_width, image_height = images[0].size

    # Calculate the size of the final grid image
    grid_width = image_width * len(x_values)
    grid_height = image_height * 1

    # Create a new blank image for the grid
    grid_image = Image.new('RGB', (grid_width, grid_height))

    # Paste images onto the grid
    for row in range(1):
        for col in range(len(x_values)):
            index = row * 4 + col
            if index < len(images):
                grid_image.paste(images[index], (col * image_width, row * image_height))

    # Save the grid image as a PNG file
    png_path = os.path.join(pdf_directory, f'image_grid_{X}.png')
    if fixed_option == "fixed": png_path = os.path.join(pdf_directory, f'image_grid_{X}_fixed.png')
    grid_image.save(png_path)

    # Clean up: Remove temporary PNG files
    for idiri in range(len(x_values)):
        png_filename = os.path.join(png_temp_directory, f"{idiri}.png")
        os.remove(png_filename)


all_dirs = [
    "./ep_10_100_K0_density_pow025/", "./ep_10_100_K0_density_pow05/", "./ep_10_100_K0_density_pow1/", "./ep_10_100_K0_density_pow15/",
    "./eC_10_100_K4_density_pow025/", "./eC_10_100_K4_density_pow05/", "./eC_10_100_K4_density_pow1/", "./eC_10_100_K4_density_pow15/",
    "./eAu_10_100_K4_density_pow025/", "./eAu_10_100_K4_density_pow05/", "./eAu_10_100_K4_density_pow1/", "./eAu_10_100_K4_density_pow15/"#,
    #"./eAu_10_100_K10_density_pow05/", "./eAu_10_100_K10_density_pow15/",
    #"./eU_10_100_K4_density_pow025/", "./eU_10_100_K4_density_pow05/", "./eU_10_100_K4_density_pow1/", "./eU_10_100_K4_density_pow15/"
]

dir = sys.argv[1]
gif_option = sys.argv[2]
fixed_option = sys.argv[3]

if gif_option == "gif":
    if dir == "grid":
        animate_grid(all_dirs, fixed_option)
    """
    else:
        animate_dir(dir, fixed_option)
        """
else:
    if dir == "all":
        for xdir in all_dirs:
            RL_row(xdir, fixed_option)
    else:
        RL_row(dir, fixed_option)
