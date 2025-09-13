// src/texture.rs
use crate::vector3::Vector3;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Texture {
    pub width: u32,
    pub height: u32,
    pub data: Vec<Vector3>, // RGB values as Vector3
}

impl Texture {
    // Create a solid color texture (fallback)
    pub fn solid_color(color: Vector3) -> Self {
        Texture {
            width: 1,
            height: 1,
            data: vec![color],
        }
    }

    // Load texture from image file
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn std::error::Error>> {
        use image::io::Reader as ImageReader;
        
        let img = ImageReader::open(path)?.decode()?;
        let rgb_img = img.to_rgb8();
        let (width, height) = rgb_img.dimensions();
        
        let mut data = Vec::with_capacity((width * height) as usize);
        
        for pixel in rgb_img.pixels() {
            let r = pixel[0] as f32 / 255.0;
            let g = pixel[1] as f32 / 255.0;
            let b = pixel[2] as f32 / 255.0;
            data.push(Vector3::new(r, g, b));
        }
        
        Ok(Texture {
            width,
            height,
            data,
        })
    }

    // Sample the texture at UV coordinates (u, v should be in [0, 1])
    pub fn sample(&self, u: f32, v: f32) -> Vector3 {
        // Wrap UV coordinates
        let u = u.fract().abs();
        let v = v.fract().abs();

        // Convert to pixel coordinates
        let x = (u * (self.width - 1) as f32) as usize;
        let y = (v * (self.height - 1) as f32) as usize;

        // Ensure we don't go out of bounds
        let x = x.min(self.width as usize - 1);
        let y = y.min(self.height as usize - 1);

        let index = y * self.width as usize + x;
        self.data.get(index).copied().unwrap_or(Vector3::new(1.0, 0.0, 1.0)) // Magenta for errors
    }
}

// Material structure
#[derive(Debug, Clone)]
pub struct Material {
    pub texture: Option<Texture>,
    pub base_color: Vector3, // Fallback color if no texture
}

impl Material {
    pub fn with_texture(texture: Texture) -> Self {
        Material {
            texture: Some(texture),
            base_color: Vector3::new(0.8, 0.8, 0.8), // Default gray
        }
    }

    pub fn solid_color(color: Vector3) -> Self {
        Material {
            texture: None,
            base_color: color,
        }
    }

    pub fn get_color_at(&self, u: f32, v: f32) -> Vector3 {
        if let Some(ref texture) = self.texture {
            texture.sample(u, v)
        } else {
            self.base_color
        }
    }
}