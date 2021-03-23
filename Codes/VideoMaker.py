#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 13:47:48 2020

@author: karthik
"""

# importing libraries 
import os 
import cv2  
  
# Checking the current directory path 
print(os.getcwd())  

  
mean_height = 0
mean_width = 0

  
  
# Video Generating function 
def generate_video(path,Iter): 
    image_folder = path # make sure to use your folder 
    video_name = os.path.join(path, "..",'Video_Iter'+str(Iter)+'.avi')
    imagelist = sorted(os.listdir(image_folder))
    imagelist = sorted(imagelist,key=len)
    imagelist = [i for i in imagelist if(i[-3:] == "png")]
    images = [img for img in imagelist] 
     
    # Array images should only consider 
    # the image files ignoring others if any 
    print(images)  
  
    frame = cv2.imread(os.path.join(image_folder, images[0])) 
  
    # setting the frame width, height width 
    # the width, height of first image 
    height, width, layers = frame.shape   
    print(frame.shape)
  
    video = cv2.VideoWriter(video_name, cv2.VideoWriter_fourcc('F', 'M', 'P', '4') , 15, (width, height))  
  
    # Appending the images to the video one by one 
    for image in images:  
        video.write(cv2.imread(os.path.join(image_folder, image)))  
      
    # Deallocating memories taken for window creation 
    cv2.destroyAllWindows()  
    video.release()  # releasing the video generated 


if __name__ =="__main__":
    print("Calling Video maker")
    
  
#Calling the generate_video function
for i in range(2,3):
    path = os.path.join("DraftIter"+str(i),"Images")
    generate_video(path,i) 