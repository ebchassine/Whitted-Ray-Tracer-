module main
import os
import math
import gfx
import rand

const (
    scene_filenames = [
        'P03_sphere_mb',
        'P03_quad_mb',        
        'P03_balls_on_plane_mb',
        'P03_ball_on_plane_texture',
        'P03_ball_on_plane_texture_mb',
        'P03_creative_artifact',
    ]
    texture_filenames = [
        'suntexture', //0
        'venustexture', //1
        'neptunetexture', //2
        'earthtexture', //3
    ]
)

type Point     = gfx.Point
type Vector    = gfx.Vector
type Direction = gfx.Direction
type Normal    = gfx.Normal
type Ray       = gfx.Ray
type Color     = gfx.Color
type Image     = gfx.Image
type Image4    = gfx.Image4
type Point2    = gfx.Point2
type Material  = gfx.Material

type Intersection = gfx.Intersection
type Surface      = gfx.Surface
type Scene        = gfx.Scene

struct Color_point {
    x int @[required]
    y int @[required]
    color Color @[required]
}
////////////////////////////////////////////////////////////////////////////////////////

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

// Interpolation function : determins the sphere's position at a given t from frame to endframe 
// Used for motion blur feature
fn interpolate(time f64, shape Surface) gfx.Frame{
    return gfx.Frame{
        o: shape.frame.o.lerp(shape.endframe.o, time)
        x: shape.frame.x
        y: shape.frame.y
        z: shape.frame.z
    }
}

fn intersect_ray_surface(surface Surface, ray Ray, time f64) Intersection {
    shape := surface.shape
    mut frame := surface.frame
    if surface.motion_blur{
        frame = interpolate(time, surface)
    }

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
        intersection_normal := frame.o.vector_to(intersection_point).normalize() // vector // use for 
        temp_frame :=  gfx.frame_oz(intersection_point, intersection_normal.as_direction()) 

        theta := math.acos( -1 * ((temp_frame.o.y - surface.frame.o.y) / surface.radius))
        phi := math.atan2(temp_frame.o.z - surface.frame.o.z, temp_frame.o.x - surface.frame.o.x)
        
        // theta := math.acos(-1 * (temp_frame.o.y )) // / 1/radius 
        // phi := math.atan2(temp_frame.o.z, temp_frame.o.x)
        u := (phi + math.pi) / (2 * math.pi) 
        v := (math.pi - theta) / math.pi
        
        texture_point := Point2{
            x: u
            y: v
        }

        return gfx.Intersection{
            frame: gfx.frame_oz(intersection_point, intersection_normal.as_direction()) 
            surface: surface
            distance: t
            texture_point: texture_point
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
        v1_temp := Vector{
            x:b.x - a.x
            y:b.y - a.y
            z:b.z - a.z
        }
        v2_temp := Vector{
            x:c.x - a.x
            y:c.y - a.y
            z:c.z - a.z
        }

        intersection_normal := v1_temp.cross(v2_temp)

        frame.o.vector_to(intersection_point).normalize() // vector
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

fn intersect_ray_scene(scene Scene, ray Ray, time f64) Intersection {
    mut closest := gfx.no_intersection  // type is Intersection

    for surface in scene.surfaces{
        mut temp_inter := gfx.no_intersection
        temp_inter = intersect_ray_surface(surface, ray, time)
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


fn map_points(intersection Intersection) Point2{
    theta := math.acos(-1 * intersection.frame.o.y) // / 1/radius 
    phi := math.atan2(intersection.frame.o.z, intersection.frame.o.x)
    u := (phi + math.pi) / (2 * math.pi) 
    v := (math.pi - theta) / math.pi
    return Point2{
        x: u
        y: v
    }
}

fn irradiance(scene Scene, ray Ray, time f64, textures []Image4) Color {
    mut accum := gfx.black
    intersection := intersect_ray_scene(scene, ray, time)
    if intersection.miss(){
        accum = scene.background_color
        return accum
    }
    mut material := intersection.surface.material

    if intersection.surface.texture_index != -1{
        texture_image := textures[intersection.surface.texture_index]
        // texture_point := map_points(intersection)
        texture_point := intersection.texture_point
        
        scaled_point := Point2{
            x: texture_point.x * 512 
            y: texture_point.y * 512 
        }
        material = gfx.Material{
            ...material
            kd: texture_image.get_color(scaled_point.as_point2i())
            
        }
    }
    if scene.ambient_color != Color{ 0.0, 0.0, 0.0 }{ //checking if there is ambient color 
        accum += scene.background_color * material.kd
    } else{
        accum = gfx.black
    }
    
    // surface.material.(kd Color, ks Color, n f64, kr Color)
    for light in scene.lights{ // implement point light from slides 
        mut l_hat := Direction{}
        mut li := Color{}
        if light.light_type == gfx.Light_Type.direction{
            l_hat = Direction{light.frame.o.x, light.frame.o.y, light.frame.o.z }
            li = light.kl // Color
        } else {
            l_hat = intersection.frame.o.direction_to(light.frame.o) // Direction
            li = light.kl.scale(1/intersection.frame.o.distance_squared_to(light.frame.o)) // Color
        }
        v := ray.d.scale(-1) 
        h := (v + l_hat.as_vector()).as_direction()
        kd := material.kd
        ks := material.ks
        normal := intersection.frame.z
        
        // Shadow Ray 
        shadow_ray := Ray{
            e: intersection.frame.o
            d: l_hat
            t_min: 0.00001
            t_max: intersection.frame.o.distance_to(light.frame.o)
        }
        shadow_intersection := intersect_ray_scene(scene, shadow_ray, time)
        if shadow_intersection.miss(){
            
            diffuse := kd //.scale(n_dot_l)
            specular := ks.scale(math.pow((math.max(0.0, normal.dot(h))), material.n))
            accum = accum.add(  ((diffuse + specular) * li).scale(math.abs(normal.dot(l_hat)))  )
        }   
    }
    // Reflection ray 
    if material.kr != gfx.black{
        v := ray.d.scale(-1) 
        normal := intersection.frame.z 
        reflect_ray := Ray{
            e: intersection.frame.o
            d: v.reflect(normal).as_direction()
        }
        reflection_c := irradiance(scene, reflect_ray, time, textures) * (material.kd)
        accum = accum.add(reflection_c)
    }
    return accum
}

// Calculates point of intersection in which the ray hits the scene
fn q_func(o Point, w f64, h f64, d f64, u f64, v f64, x Direction, y Direction, z Direction) Point{
    return o.add(x.scale((u-0.5) * w) + y.scale((v-0.5)*h) + z.scale(-d))
}

fn raytrace_helper(scene Scene, textures []Image4, y int) []Color_point{
    w := scene.camera.sensor.resolution.width
    h := scene.camera.sensor.resolution.height
    mut column := []Color_point{}
    for x in 0 .. w{
        u := (f64(x)/f64(w))
        v := (f64(y)/f64(h))
        q := q_func(scene.camera.frame.o, scene.camera.sensor.size.width, scene.camera.sensor.size.height, scene.camera.sensor.distance, u, v, scene.camera.frame.x, scene.camera.frame.y, scene.camera.frame.z)
        
        curr_ray := Ray(scene.camera.frame.o.ray_through(q)) 
        time := rand.f64_in_range(0.0, 1.0) or {0}  // Create random time here 

        curr_color := irradiance(scene, curr_ray, time, textures) 

        column << Color_point{
            x: x
            y: y
            color: curr_color
        }
    }
    return column  
}

fn raytrace_helper_aa(scene Scene, textures []Image4, y int, samples int) []Color_point{
    w := scene.camera.sensor.resolution.width
    h := scene.camera.sensor.resolution.height
    mut column := []Color_point{}
    for x in  0 .. w{
        mut accum := Color{}
        for y_sample in 0 .. samples{
            for x_sample in 0 .. samples{
                u := (f64(x) + (f64(x_sample) + 0.5) / f64(samples)) / f64(w)
                v := (f64(y) + (f64(y_sample) + 0.5) / f64(samples)) / f64(h)
                q := q_func(scene.camera.frame.o, scene.camera.sensor.size.width, scene.camera.sensor.size.height, scene.camera.sensor.distance, u, v, scene.camera.frame.x, scene.camera.frame.y, scene.camera.frame.z)
                curr_ray := Ray(scene.camera.frame.o.ray_through(q))
                time := rand.f64_in_range(0.0, 1.0) or {0}
                temp_c := irradiance(scene, curr_ray, time, textures)
                accum = accum.add(temp_c)
            }
        }
        final_c := accum.scale(1.0 / f64(samples * samples))
        column << Color_point{
            x: x
            y: h - y
            color: final_c
        }
    }
    return column
}

// Computes image of scene using basic Whitted raytracer.
fn raytrace(scene Scene, textures []Image4) Image { // pass list of Image types as textures 
    // w := scene.camera.sensor.resolution.width
    h := scene.camera.sensor.resolution.height
    mut image := gfx.Image.new(scene.camera.sensor.resolution)
    mut threads := []thread []Color_point{}

    if scene.camera.sensor.samples == 1{
        for y in 0 .. h{
            // in this case, each thread is an array filled with Color_points, with each representing a column's colors
            threads << spawn raytrace_helper(scene, textures, y)  
            // for x in  0 .. w{ This is where the original loop occurred 
            // }
        }
        println("GOT HERE")
        res := threads.wait()
        // println('Threads finished: ${res}')
        for column in res{
            for pixel in column{
                x := pixel.x
                y := pixel.y
                color := pixel.color
                image.set_xy(x, y, color)
            }
        }
    } else {
        samples := scene.camera.sensor.samples
        for y in 0 .. h{
            threads << spawn raytrace_helper_aa(scene, textures, y, samples)
        }
        res := threads.wait()
        for column in res{
            for pixel in column{
                x := pixel.x
                y := pixel.y
                color := pixel.color
                image.set_xy(x, y, color)
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
    // Loading all textures from texture folder into array to be passed into raytracer -> irradiance
    mut texture_array := []Image4{}
    for filename in texture_filenames{
        println('Rendering texture ${filename}...')
        texture_array << gfx.load_png('textures/${filename}.png')
    }

    for filename in scene_filenames {
        println('Rendering ${filename}...')
        scene := gfx.scene_from_file('P03_Final/scenes/${filename}.json')!
        image := raytrace(scene, texture_array)
        image.save_png('output/${filename}.png')
    }
    println('Done!')
}