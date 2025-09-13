use crate::vector3::Vector3;

#[derive(Debug, Clone, Copy)]
pub struct DirectionalLight {
    pub direction: Vector3,  // Direction TO the light (normalized)
    pub color: Vector3,
    pub intensity: f32,
}

impl DirectionalLight {
    pub fn new(direction: Vector3, color: Vector3, intensity: f32) -> Self {
        DirectionalLight {
            direction: (-direction).normalize(), // We want direction TO light
            color,
            intensity,
        }
    }

    pub fn calculate_diffuse(&self, surface_normal: Vector3, surface_color: Vector3) -> Vector3 {
        // Lambertian diffuse shading
        let cos_theta = surface_normal.dot(self.direction).max(0.0);
        let diffuse_strength = cos_theta * self.intensity;
        
        Vector3::new(
            surface_color.x * self.color.x * diffuse_strength,
            surface_color.y * self.color.y * diffuse_strength,
            surface_color.z * self.color.z * diffuse_strength,
        )
    }
}

#[derive(Debug, Clone, Copy)]
pub struct AmbientLight {
    pub color: Vector3,
    pub intensity: f32,
}

impl AmbientLight {
    pub fn new(color: Vector3, intensity: f32) -> Self {
        AmbientLight { color, intensity }
    }

    pub fn calculate_ambient(&self, surface_color: Vector3) -> Vector3 {
        Vector3::new(
            surface_color.x * self.color.x * self.intensity,
            surface_color.y * self.color.y * self.intensity,
            surface_color.z * self.color.z * self.intensity,
        )
    }
}