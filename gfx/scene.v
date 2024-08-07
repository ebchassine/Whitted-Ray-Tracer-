module gfx

import os
import json


///////////////////////////////////////////////////////////
// structs used to define a scene
//
// NOTE: if any post-load processing needs to happen on
//       structs after loading from JSON file, add
//       processing code to `update_after_load` below!
//
//       ex: Frames can have Look At parameters (eye, target, up)
//           specified instead of usual Frame params (o, x, y, z).

pub struct Scene {
pub:
    camera           Camera
    background_color Color       = Color{ 0.2, 0.2, 0.2 }
    ambient_color    Color       = Color{ 0.2, 0.2, 0.2 }
    lights           [] Light    = [ Light{} ]
    surfaces         [] Surface  = [ Surface{} ]
    accel            [] Accel    = []
}

pub struct Accel {
pub:
    radius f64  = 1.0
    frame Frame = Frame{
        o: Point{ 0.0, 0.0, 1.0 }
        x: direction_x
        y: direction_y
        z: direction_z
    }
    meshes [] TriangleMesh = []
}

pub struct TriangleMesh {
pub:
    triangles [] Triangle = []
}

pub struct Triangle {
pub:
    p0  Point  = Point{  0.0 0.0 0.0 }
    kd0 Color  = gfx.white
    uv0 Point2 = Point2{ 0.0 0.0 }

    p1  Point  = Point{  1.0 0.0 0.0 }
    kd1 Color  = gfx.white
    uv1 Point2 = Point2{ 1.0 0.0 }

    p2  Point  = Point{  0.0 1.0 0.0 }
    kd2 Color  = gfx.white
    uv2 Point2 = Point2{ 0.0 1.0 }
}

pub struct Camera {
pub:
    sensor Sensor = Sensor{}
    frame  Frame  = Frame{
        o: Point{ 0.0, 0.0, 1.0 }
        x: direction_x
        y: direction_y
        z: direction_z
    }
}

pub struct Sensor {
pub:
    size       Size2  = Size2{  1.0, 1.0 }
    resolution Size2i = Size2i{ 512, 512 }
    distance   f64    = 1.0
    samples    int    = 1
}
pub enum Light_Type{
    point
    direction
}
pub struct Light {
pub:
    kl     Color = white
    frame  Frame = Frame{
        o: Point{ 0.0, 0.0, 5.0 }
        x: direction_x
        y: direction_y
        z: direction_z
    }
    light_type Light_Type = Light_Type.point
}

pub enum Shape {
    sphere
    quad
    triangle
    circle
}

pub struct Surface {
pub:
    shape    Shape = Shape.sphere
    radius   f64   = 1.0
    frame    Frame
    material Material
}

pub struct Material {
pub:
    kd Color = white
    ks Color = black
    n  f64   = 10.0
    kr Color = black
}


// function to update scene after loading from JSON
fn update_after_load(scene Scene) Scene {
    // NOTE: attempted a more generic version using reflection (see bottom)
    return Scene{
        ...scene                                            // copy everything from scene with some overrides (below)
        camera: Camera{
            ...scene.camera                                 // copy everything from camera with override
            frame: scene.camera.frame.as_frame()            // update frame for camera if look at is specified
        }
        lights: scene.lights.map(fn (light Light) Light {   // for each light...
            return Light{
                ...light                                    // copy everything from camera with override
                frame: light.frame.as_frame()               // update frame for light if look at is specified
            }
        })

    }
}



///////////////////////////////////////////////////////////
// scene importing and exporting functions

pub fn scene_from_file(path string) !Scene {
    // load and decode scene from JSON file
    data := os.read_file(path)!
    scene := json.decode(Scene, data)!
    return update_after_load(scene)
}

pub fn (scene Scene) to_json() string {
    return json.encode(scene)
}



///////////////////////////////////////////////////////////
// convenience getters

pub fn (light Light) o() Point {
    return light.frame.o
}

pub fn (surface Surface) o() Point {
    return surface.frame.o
}



///////////////////////////////////////////////////////////
// attempting a more generic updated function that uses reflection

// // commented out following code due to the wrong struct data is seen on the lines marked with <---
// type SceneType = Scene | Camera | Sensor | Light | Surface | Material
// fn update_frame_from_load(f Frame) Frame {
//     println('UPDATING FRAME FROM LOAD')
//     println(f)                                          // <---  this has INCORRECT data!!
//     if f.eye != f.target && f.up != direction_no {
//         // look at was specified!
//         ret := frame_lookat(f.eye, f.target, f.up)
//         println(ret)
//         return ret
//     }
//     return f
// }
// fn update_array_from_load[T](a []T) []T {
//     mut result := []T{}
//     for i in a {
//         result << update_after_load(i)
//     }
//     return result
// }
// fn update_after_load[T](o T) T {
//     mut result := T{}
//     $for field in T.fields {
//         $if field.typ is $array {
//             result.$(field.name) = update_array_from_load(o.$(field.name))
//         } $else $if field.typ is $enum {
//             result.$(field.name) = o.$(field.name)
//         } $else $if field.typ is Frame {
//             println(o.$(field.name))                                        // <---  this has INCORRECT data!!
//             result.$(field.name) = update_frame_from_load(o.$(field.name))
//         } $else $if field.typ is SceneType {
//             result.$(field.name) = update_after_load(o.$(field.name))
//         } $else {
//             result.$(field.name) = o.$(field.name)
//         }
//     }
//     return result
// }





