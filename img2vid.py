#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 20:05:49 2023

@author: tadewuyi
"""

import os
import cv2
from PIL import Image


def img2vid(path: str, filename: str, ext: str = '.png'):
  '''
  
  '''
  images = [img for img in os.listdir(path) if img.endswith(ext)]
  
  images.sort()
  output_path = "../"
  video_filename = f'{filename}.mp4'
  video_path = os.path.join(output_path, video_filename)
#  codec = cv2.VideoWriter_fourcc(*'mp4v')
#  
#  # Open the video writer
#  frame = cv2.imread(os.path.join(path, images[0]))
#  height, width, layers = frame.shape
#  video = cv2.VideoWriter(video_path, codec, 1, (width, height))
#  
#  # Write images to the video
#  for image in images:
#      video.write(Image.open(os.path.join(path, image)))
#  
#  # Close the video writer
#  cv2.destroyAllWindows()
#  video.release()
  
  
  
  
  first_image = Image.open(os.path.join(path, images[0]))
  width, height = first_image.size
  image_objects = []
  
  
  for img in images:
    image_objects.append([Image.open(os.path.join(path, img))])
    
#  image_objects = [Image.open(os.path.join(path, img)) for img in images]
  
  image_objects[0].save(
      video_path,
      save_all=True,
      append_images=image_objects[1:],
      duration=5,  # Specify the duration between frames in milliseconds
      loop=0)  # Specify the number of loops (0 means infinite loop)

  
  print("Video creation completed!")








































