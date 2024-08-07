module main
import os
import math
import gfx

////////////////////////////////////////////////////////////////////////////////////////
// Comment out lines in array below to prevent re-rendering every scene.
// If you create a new scene file, add it to the list below.
//
// NOTE: **BEFORE** you submit your solution, uncomment all lines, so
//       your code will render all the scenes!

const (
    scene_filenames = [
        'P02_00_sphere',
        'P02_01_sphere_ambient',
        'P02_02_sphere_room',
        'P02_03_quad',
        'P02_04_quad_room',
        'P02_05_ball_on_plane',
        'P02_06_balls_on_plane',
        'P02_07_reflections',
        'P02_08_antialiased',
        'P02_09_triangle',
        'P02_09_planar_circle',
        'P02_09_cube',
        'P02_09_ball_on_plane_light',
        'P02_10_quad_room_artifact',
    ]
)

////////////////////////////////////////////////////////////////////////////////////////
// module aliasing to make code a little easier to read
// ex: replacing `gfx.Scene` with just `Scene`

type Point     = gfx.Point
type Vector    = gfx.Vector
type Direction = gfx.Direction
type Normal    = gfx.Normal
type Ray       = gfx.Ray
type Color     = gfx.Color
type Image     = gfx.Image

type Intersection = gfx.Intersection
type Surface      = gfx.Surface
type Scene        = gfx.Scene

////////////////////////////////////////////////////////////////////////////////////////
// functions to implement

// Helper functions for point arithemtic, I made these during triangle-ray intersect which is why they are only used there
fn dot_prod(a Point, b Point) f64{
    return a.x * b.x + a.y * b.y + a.z * b.z 
}
fn scale_point(a Point, b f64) Point{
    return Point{
        x: a.x * b
        y: a.y * b
        z: a.z * b
    }
}
fn cross_point(a Point, b Point ) Point {
    return Point{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    }
}

fn intersect_ray_surface(surface Surface, ray Ray) Intersection {
    shape := surface.shape
    frame := surface.frame

    if shape == gfx.Shape.sphere{
        sphere_center := frame.o
        r := surface.radius // f64
        e := ray.e // Point 
        d := ray.d

        a := ray.d.dot(ray.d)
        e_sub_c := (Point{ e.x - sphere_center.x, e.y - sphere_center.y, e.z - sphere_center.z } )
        b := 2 * (d.dot(e_sub_c))
        c :=  e_sub_c.x * e_sub_c.x + e_sub_c.y * e_sub_c.y + e_sub_c.z * e_sub_c.z - r * r
        discriminant := b * b - (4 * a * c)

        if discriminant  < 0 {
            return gfx.no_intersection
        }

        t_1 := ((-b - math.sqrt(discriminant)) / (2 * a))
        t_2 := ((-b + math.sqrt(discriminant)) / (2 * a))
        mut t := 0.0
        if t_1 > ray.t_min && t_1 < ray.t_max{
            t = t_1
        } else if  t_2 > ray.t_min && t_2 < ray.t_max{
            t = t_2
        } else {
            return gfx.no_intersection
        }
        
        intersection_point := ray.at(t) // point 
        intersection_normal := frame.o.vector_to(intersection_point).normalize() // vector
        return gfx.Intersection{
            frame: gfx.frame_oz(intersection_point, intersection_normal.as_direction()) //check?
            surface: surface
            distance: t
        }
    } else if shape == gfx.Shape.quad{
        center := frame.o // point 
        normal := frame.z // direction
        e := ray.e // point 
        d_hat := ray.d // direction

        c_sub_e := Point{ center.x - e.x, center.y - e.y, center.z - e.z }
        numer := normal.dot(c_sub_e)
        denom := d_hat.dot(normal)

        t := numer / denom
        if t < ray.t_min || t > ray.t_max {
            return gfx.no_intersection
        }
        intersection_point := ray.at(t)
        intersection_normal := normal // or surface.frame.z
        loc_point := Direction{ intersection_point.x - center.x, intersection_point.y - center.y, intersection_point.z - center.z}
        if loc_point.linf_norm() > surface.radius{
            return gfx.no_intersection
        }
        return gfx.Intersection{    
            frame: gfx.frame_oz(intersection_point, intersection_normal) //check?
            surface: surface
            distance: t
        }
    } else if shape == gfx.Shape.triangle{
        a := Point{
            x:surface.frame.x.x
            y:surface.frame.x.y
            z:surface.frame.x.z 
        }
        b := Point {
            x:surface.frame.y.x
            y:surface.frame.y.y
            z:surface.frame.y.z 
        }
        c := Point {
            x:surface.frame.z.x
            y:surface.frame.z.y
            z:surface.frame.z.z 
        }
        e_prime := Point{ray.e.x - c.x, ray.e.y - c.x, ray.e.z - c.z} // Point 
        a_prime := Point{a.x - c.x, a.y - c.x, a.z - c.z} // Point 
        b_prime := Point{b.x - c.x, b.y - c.x, b.z - c.z} // Point 
        d := Point{ 
            x:ray.d.x
            y:ray.d.y
            z:ray.d.z
        }
        e_cross_a := cross_point(e_prime, a)
        d_cross_b := cross_point(d, b_prime)
        t_num := dot_prod(e_cross_a, b)
        t_denom := dot_prod(d_cross_b, a)
        t := t_num / t_denom

        alpha_num := dot_prod(d_cross_b, e_prime)
        alpha_denom := dot_prod(d_cross_b, a_prime)
        alpha := alpha_num / alpha_denom

        beta_num := dot_prod(e_cross_a, d)
        beta_denom := dot_prod(d_cross_b, a_prime)
        beta := beta_num / beta_denom

        if t < ray.t_min || t > ray.t_max || alpha < 0 || beta < 0 || alpha + beta > 1  {
            return gfx.no_intersection
        }

        intersection_point := ray.at(t)
        intersection_normal := frame.o.vector_to(intersection_point).normalize() // vector
        return gfx.Intersection{
            frame: gfx.frame_oz(intersection_point, intersection_normal.as_direction()) //check?
            surface: surface
            distance: t
        }
    
    } else if shape == gfx.Shape.circle{
        center := frame.o // point 
        normal := frame.z // direction
        e := ray.e // point 
        d_hat := ray.d // direction

        c_sub_e := Point{ center.x - e.x, center.y - e.y, center.z - e.z }
        numer := normal.dot(c_sub_e)
        denom := d_hat.dot(normal)

        t := numer / denom
        if t < ray.t_min || t > ray.t_max {
            return gfx.no_intersection
        }
        intersection_point := ray.at(t)
        intersection_normal := normal // or surface.frame.z
        loc_point := Direction{ intersection_point.x - center.x, intersection_point.y - center.y, intersection_point.z - center.z}
        if loc_point.l2_norm() > surface.radius{
            return gfx.no_intersection
        }
        return gfx.Intersection{    
            frame: gfx.frame_oz(intersection_point, intersection_normal) //check?
            surface: surface
            distance: t
        }
    }

    return gfx.no_intersection
}


fn intersect_ray_scene(scene Scene, ray Ray) Intersection {
    mut closest := gfx.no_intersection  // type is Intersection

    for surface in scene.surfaces{
        temp_inter := intersect_ray_surface(surface, ray)

        if temp_inter.miss(){
            continue 
        // } if temp_inter.distance() > closest.distance(){
            // continue
        } else if closest.miss() || temp_inter.distance <= closest.distance + 1e-6{
            closest = temp_inter
        } 
    }
    return closest  
}

fn irradiance(scene Scene, ray Ray) Color {
    mut accum := gfx.black

    intersection := intersect_ray_scene(scene, ray)
    if scene.ambient_color != Color{ 0.0, 0.0, 0.0 }{ //checking if there is ambient color 
        accum = scene.background_color * intersection.surface.material.kd
    } else{
        accum = gfx.black
    }
    
    if intersection.miss(){
        accum = scene.background_color
        return accum
    } 
    
    // surface.material.(kd Color, ks Color, n f64, kr Color)
    for light in scene.lights{ // implement point light from slides 

        mut l_hat := Direction{}
        mut li := Color{}
        if light.light_type == gfx.Light_Type.direction{
            l_hat = Direction{light.frame.o.x, light.frame.o.y, light.frame.o.z }
            // li = light.kl
            li = light.kl.scale(1/intersection.frame.o.distance_squared_to(light.frame.o)) // Color
        } else {
            l_hat = intersection.frame.o.direction_to(light.frame.o) // Direction
            li = light.kl.scale(1/intersection.frame.o.distance_squared_to(light.frame.o)) // Color
        }
        // l_hat := intersection.frame.o.direction_to(light.frame.o) // Direction
        // li := light.kl.scale(1/intersection.frame.o.distance_squared_to(light.frame.o)) // Color
        v := ray.d.scale(-1) 
        h := (v + l_hat.as_vector()).as_direction()
        kd := intersection.surface.material.kd
        ks := intersection.surface.material.ks
        normal := intersection.frame.z
        
        // Shadow Ray 
        shadow_ray := Ray{
            e: intersection.frame.o
            d: l_hat
            t_min: 0.00001
            t_max: intersection.frame.o.distance_to(light.frame.o)
        }
        shadow_intersection := intersect_ray_scene(scene, shadow_ray)
        if shadow_intersection.miss(){
            diffuse := kd //.scale(n_dot_l)
            specular := ks.scale(math.pow((math.max(0.0, normal.dot(h))), intersection.surface.material.n))
            accum = accum.add(  ((diffuse + specular) * li).scale(math.abs(normal.dot(l_hat)))  )
        }   
    }
    // Reflection ray 
    if intersection.surface.material.kr != gfx.black{
        v := ray.d.scale(-1) 
        normal := intersection.frame.z
        reflect_ray := Ray{
            e: intersection.frame.o
            d: v.reflect(normal).as_direction()
        }
        reflection_c := irradiance(scene, reflect_ray) * (intersection.surface.material.kd)
        accum = accum.add(reflection_c)
    }

    return accum
}

fn q_func(o Point, w f64, h f64, d f64, u f64, v f64, x Direction, y Direction, z Direction) Point{
    return o.add(x.scale((u-0.5) * w) + y.scale((v-0.5)*h) + z.scale(-d))
}
// Computes image of scene using basic Whitted raytracer.
fn raytrace(scene Scene) Image {
    // println(scene)
    w := scene.camera.sensor.resolution.width
    h := scene.camera.sensor.resolution.height
    mut image := gfx.Image.new(scene.camera.sensor.resolution)
    if scene.camera.sensor.samples == 1{ //A_a line 
        for y in 0 .. h{
            for x in  0 .. w{
                u := (f64(x)/f64(w))
                v := (f64(y)/f64(h))
                q := q_func(scene.camera.frame.o, scene.camera.sensor.size.width, scene.camera.sensor.size.height, scene.camera.sensor.distance, u, v, scene.camera.frame.x, scene.camera.frame.y, scene.camera.frame.z)
                curr_ray := Ray(scene.camera.frame.o.ray_through(q))
                curr_color := irradiance(scene, curr_ray)
                image.set_xy(x, h - y, curr_color)
            }
        }
    } else {
        samples := scene.camera.sensor.samples
        for y in 0 .. h{
            for x in  0 .. w{
                mut accum := Color{}
                for y_sample in 0 .. samples{
                    for x_sample in 0 .. samples{
                        u := (f64(x) + (f64(x_sample) + 0.5) / f64(samples)) / f64(w)
                        v := (f64(y) + (f64(y_sample) + 0.5) / f64(samples)) / f64(h)
                        q := q_func(scene.camera.frame.o, scene.camera.sensor.size.width, scene.camera.sensor.size.height, scene.camera.sensor.distance, u, v, scene.camera.frame.x, scene.camera.frame.y, scene.camera.frame.z)
                        curr_ray := Ray(scene.camera.frame.o.ray_through(q))
                        temp_c := irradiance(scene, curr_ray)
                        accum = accum.add(temp_c)
                    }
                }
                final_c := accum.scale(1.0 / f64(samples * samples))
                image.set_xy(x, h - y, final_c)
            }
        }
    }
    return image
}

fn main() {
    // Make sure images folder exists, because this is where all generated images will be saved
    if !os.exists('output') {
        os.mkdir('output') or { panic(err) }
    }

    for filename in scene_filenames {
        println('Rendering ${filename}...')
        scene := gfx.scene_from_file('P02_Raytrace/scenes/${filename}.json')!
        image := raytrace(scene)
        image.save_png('output/${filename}.png')
    }

    println('Done!')
}