// src/raytracer.rs
use crate::vector3::Vector3;
use crate::ray::{Ray, HitRecord};
use crate::camera::Camera;
use crate::cube::{Cube, Plane};
use crate::light::{DirectionalLight, AmbientLight};
use crate::texture::{Texture, Material};

pub struct Scene {
    pub cube: Cube,
    pub ground: Plane,
    pub directional_light: DirectionalLight,
    pub ambient_light: AmbientLight,
    pub background_color: Vector3,
}

impl Scene {
    pub fn new() -> Self {
        // Try to load the cube texture from image
        let cube_material = match Texture::from_file("cube_texture.jpg")
        {
            Ok(texture) => {
                println!("Texture loaded successfully!");
                Material::with_texture(texture)
            },
            Err(e) => {
                println!("Could not load texture image: {}", e);
                println!("Make sure you have a file named 'cube_texture.jpg' or 'texture.jpg' in your project root");
                println!("Supported formats: .jpg, .png, .bmp, .tiff");
                println!("Using solid red color as fallback");
                Material::solid_color(Vector3::new(0.8, 0.2, 0.2))
            }
        };

        Scene {
            // Cube with image texture (or solid color fallback)
            cube: Cube::new(
                Vector3::new(-0.5, 0.0, -0.5),
                Vector3::new(0.5, 1.0, 0.5),
                cube_material,
            ),
            // Simple colored ground (no texture)
            ground: Plane::new(
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(0.0, 1.0, 0.0),
                Vector3::new(0.65, 0.7, 0.45), // Olive green
            ),
            // Directional light from upper left
            directional_light: DirectionalLight::new(
                Vector3::new(1.0, -1.0, -1.0),
                Vector3::new(1.0, 1.0, 1.0),
                0.8,
            ),
            // Ambient light for softer shadows
            ambient_light: AmbientLight::new(
                Vector3::new(1.0, 1.0, 1.0),
                0.2,
            ),
            // Background color
            background_color: Vector3::new(0.6, 0.6, 0.3),
        }
    }

    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        let mut closest_hit: Option<HitRecord> = None;
        let mut closest_t = t_max;

        // Check cube intersection
        if let Some(hit) = self.cube.hit(ray, t_min, closest_t) {
            closest_t = hit.t;
            closest_hit = Some(hit);
        }

        // Check ground intersection
        if let Some(hit) = self.ground.hit(ray, t_min, closest_t) {
            closest_t = hit.t;
            closest_hit = Some(hit);
        }

        closest_hit
    }

    fn is_in_shadow(&self, point: Vector3, light_direction: Vector3) -> bool {
        // Create shadow ray from hit point towards light
        let shadow_ray = Ray::new(point + light_direction * 0.001, light_direction);
        
        // Check if shadow ray hits any object
        if let Some(_) = self.cube.hit(&shadow_ray, 0.001, f32::INFINITY) {
            return true;
        }

        false
    }

    fn calculate_lighting(&self, hit: &HitRecord) -> Vector3 {
        // Start with ambient lighting
        let ambient = self.ambient_light.calculate_ambient(hit.material_color);
        
        // Check if point is in shadow
        let in_shadow = self.is_in_shadow(hit.point, self.directional_light.direction);
        
        let diffuse = if in_shadow {
            // If in shadow, reduce diffuse lighting significantly
            self.directional_light.calculate_diffuse(hit.normal, hit.material_color) * 0.1
        } else {
            // Full diffuse lighting
            self.directional_light.calculate_diffuse(hit.normal, hit.material_color)
        };

        ambient + diffuse
    }
}

pub struct RayTracer {
    pub width: u32,
    pub height: u32,
    pub camera: Camera,
    pub scene: Scene,
}

impl RayTracer {
    pub fn new(width: u32, height: u32) -> Self {
        let aspect_ratio = width as f32 / height as f32;
        
        // Position camera to view the cube from an angle
        let camera = Camera::new(
            Vector3::new(3.0, 1.5, 2.5),
            Vector3::new(0.0, 0.5, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            45.0,
            aspect_ratio,
        );

        RayTracer {
            width,
            height,
            camera,
            scene: Scene::new(),
        }
    }

    pub fn trace_ray(&self, ray: &Ray) -> Vector3 {
        if let Some(hit) = self.scene.hit(ray, 0.001, f32::INFINITY) {
            self.scene.calculate_lighting(&hit)
        } else {
            self.scene.background_color
        }
    }

    pub fn render(&self) -> Vec<u8> {
        let mut pixels = vec![0u8; (self.width * self.height * 4) as usize];
        
        for y in 0..self.height {
            for x in 0..self.width {
                let u = x as f32 / (self.width - 1) as f32;
                let v = (self.height - 1 - y) as f32 / (self.height - 1) as f32; // Flip Y
                
                let ray = self.camera.get_ray(u, v);
                let color = self.trace_ray(&ray);
                
                // Clamp and convert to 0-255 range
                let r = (color.x.min(1.0).max(0.0) * 255.0) as u8;
                let g = (color.y.min(1.0).max(0.0) * 255.0) as u8;
                let b = (color.z.min(1.0).max(0.0) * 255.0) as u8;
                
                let pixel_index = ((y * self.width + x) * 4) as usize;
                pixels[pixel_index] = r;
                pixels[pixel_index + 1] = g;
                pixels[pixel_index + 2] = b;
                pixels[pixel_index + 3] = 255; // Alpha
            }
            
            // Print progress
            if y % (self.height / 10) == 0 {
                println!("Rendering: {}%", (y * 100) / self.height);
            }
        }
        
        println!("Rendering: 100%");
        pixels
    }
}