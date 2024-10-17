### Whitted Ray Tracer 


This project uses light-ray intersection calculations based on the .json file that specifies shapes (coordinates, start and endframe, and surface/texture qualities), camera (coordinates, orientation, resolution), and light sources (type, coordinate, orientation).  By shooting light-rays based on how the light source is configured, the rays may intersect (or not) shapes or reflective surfaces.  At this point the ray and surface it intersects will determine the color of the pixel painted on that x,y on the image.

The project simulates a simple ray tracer in VLang, a low level language with garbage collection, optimized for graphics processing,
and also optimized syntax for some ease of use. 

The project itself renders scenes from .json files and runs it through a pipeline in P02_Raytrace  
