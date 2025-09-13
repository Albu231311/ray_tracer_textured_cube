use crate::vector3::Vector3;

#[derive(Debug, Clone, Copy)]
pub struct Ray {
    pub origin: Vector3,
    pub direction: Vector3,
}

impl Ray {
    pub fn new(origin: Vector3, direction: Vector3) -> Self {
        Ray { 
            origin, 
            direction: direction.normalize() 
        }
    }

    pub fn at(&self, t: f32) -> Vector3 {
        self.origin + self.direction * t
    }
}

#[derive(Debug, Clone, Copy)]
pub struct HitRecord {
    pub point: Vector3,
    pub normal: Vector3,
    pub t: f32,
    pub front_face: bool,
    pub material_color: Vector3,
}

impl HitRecord {
    pub fn new(point: Vector3, normal: Vector3, t: f32, material_color: Vector3) -> Self {
        HitRecord {
            point,
            normal,
            t,
            front_face: true,
            material_color,
        }
    }

    pub fn set_face_normal(&mut self, ray: &Ray, outward_normal: Vector3) {
        self.front_face = ray.direction.dot(outward_normal) < 0.0;
        self.normal = if self.front_face { 
            outward_normal 
        } else { 
            -outward_normal 
        };
    }
}