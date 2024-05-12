# -*- coding: utf-8 -*-
"""
Created on Wed May  1 00:16:31 2024

@author: may7e
"""
'''
This script basically plots collages of selected ENA temperature maps along with their respective network 
diagrams. This is an ad hoc script and doesn't really do thing automatically. You'd need to move the images you
want to put into a collage into a separate folder and pass the path to the folder into the function. 
The part of the function that says "my_order" is to arrange the collage into a specific way. Comment that section 
out if it isn't required'. This also doesn't save the plot, just displays it. 
'''


import os
import matplotlib.pyplot as plt
from PIL import Image


def img_collage(folder_path):
    # Get list of image file names in the folder
    image_files = [file for file in os.listdir(folder_path) if file.endswith('.jpg') or file.endswith('.png')]

    # Limit to first 16 images if more than 16 are present
    image_files = image_files[:16]
    # my_order = [8,9,10,11,0,1,2,3,4,5,6,7,12,13,14,15]
    
    # image_files = [image_files[i] for i in my_order]

    # Define the grid size
    num_rows = 4
    num_cols = 4

    # Create subplots for the images
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(12, 12))

    # Iterate over images and plot in the grid
    for i, image_file in enumerate(image_files):
        # Open the image file using PIL
        image_path = os.path.join(folder_path, image_file)
        image = Image.open(image_path)

        # Plot the image on the corresponding subplot
        row = i // num_cols
        col = i % num_cols
        axs[row, col].imshow(image)
        axs[row, col].axis('off')  # Turn off axes

    plt.tight_layout()
    plt.show()