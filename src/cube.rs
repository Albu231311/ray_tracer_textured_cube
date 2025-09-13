// src/cube.rs
use crate::vector3::Vector3;
use crate::ray::{Ray, HitRecord};
use crate::texture::Material;

#[derive(Debug, Clone)]
pub struct Cube {
    pub min: Vector3,
    pub max: Vector3,
    pub material: Material,
}

impl Cube {
    pub fn new(min: Vector3, max: Vector3, material: Material) -> Self {
        Cube { min, max, material }
    }

    // Calculate UV coordinates for each face of the cube
    fn calculate_uv(&self, hit_point: Vector3, face: u32) -> (f32, f32) {
        let size = self.max - self.min;
        
        match face {
            0 => { // -X face (left)
                let u = (hit_point.z - self.min.z) / size.z;
                let v = (hit_point.y - self.min.y) / size.y;
                (u, 1.0 - v)
            },
            1 => { // +X face (right)
                let u = (self.max.z - hit_point.z) / size.z;
                let v = (hit_point.y - self.min.y) / size.y;
                (u, 1.0 - v)
            },
            2 => { // -Y face (bottom)
                let u = (hit_point.x - self.min.x) / size.x;
                let v = (hit_point.z - self.min.z) / size.z;
                (u, v)
            },
            3 => { // +Y face (top)
                let u = (hit_point.x - self.min.x) / size.x;
                let v = (self.max.z - hit_point.z) / size.z;
                (u, v)
            },
            4 => { // -Z face (back)
                let u = (self.max.x - hit_point.x) / size.x;
                let v = (hit_point.y - self.min.y) / size.y;
                (u, 1.0 - v)
            },
            5 => { // +Z face (front)
                let u = (hit_point.x - self.min.x) / size.x;
                let v = (hit_point.y - self.min.y) / size.y;
                (u, 1.0 - v)
            },
            _ => (0.0, 0.0),
        }
    }

    pub fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        // AABB (Axis-Aligned Bounding Box) intersection
        let mut t_near = t_min;
        let mut t_far = t_max;
        let mut normal = Vector3::zero();
        let mut hit_face = 0;

        // Check intersection with each pair of planes
        for axis in 0..3 {
            let (min_val, max_val, dir_val) = match axis {
                0 => (self.min.x, self.max.x, ray.direction.x),
                1 => (self.min.y, self.max.y, ray.direction.y),
                2 => (self.min.z, self.max.z, ray.direction.z),
                _ => unreachable!(),
            };

            let origin_val = match axis {
                0 => ray.origin.x,
                1 => ray.origin.y,
                2 => ray.origin.z,
                _ => unreachable!(),
            };

            if dir_val.abs() < 1e-8 {
                // Ray is parallel to the planes
                if origin_val < min_val || origin_val > max_val {
                    return None;
                }
            } else {
                let inv_dir = 1.0 / dir_val;
                let mut t0 = (min_val - origin_val) * inv_dir;
                let mut t1 = (max_val - origin_val) * inv_dir;
                let mut face0 = axis * 2;
                let mut face1 = axis * 2 + 1;

                if inv_dir < 0.0 {
                    std::mem::swap(&mut t0, &mut t1);
                    std::mem::swap(&mut face0, &mut face1);
                }

                if t0 > t_near {
                    t_near = t0;
                    hit_face = face0;
                }

                if t1 < t_far {
                    t_far = t1;
                }

                if t_near > t_far {
                    return None;
                }
            }
        }

        if t_near >= t_min && t_near <= t_max {
            let hit_point = ray.at(t_near);
            
            // Calculate normal based on which face was hit
            normal = match hit_face {
                0 => Vector3::new(-1.0, 0.0, 0.0), // -X face
                1 => Vector3::new(1.0, 0.0, 0.0),  // +X face
                2 => Vector3::new(0.0, -1.0, 0.0), // -Y face
                3 => Vector3::new(0.0, 1.0, 0.0),  // +Y face
                4 => Vector3::new(0.0, 0.0, -1.0), // -Z face
                5 => Vector3::new(0.0, 0.0, 1.0),  // +Z face
                _ => Vector3::new(0.0, 1.0, 0.0),
            };

            // Calculate UV coordinates for this face
            let (u, v) = self.calculate_uv(hit_point, hit_face);
            
            // Get the color from the material at this UV coordinate
            let material_color = self.material.get_color_at(u, v);
            
            let mut hit_record = HitRecord::new(hit_point, normal, t_near, material_color);
            hit_record.set_face_normal(ray, normal);
            
            Some(hit_record)
        } else {
            None
        }
    }
}

// Simple ground plane (sin textura)
#[derive(Debug, Clone, Copy)]
pub struct Plane {
    pub point: Vector3,
    pub normal: Vector3,
    pub color: Vector3,
}

impl Plane {
    pub fn new(point: Vector3, normal: Vector3, color: Vector3) -> Self {
        Plane {
            point,
            normal: normal.normalize(),
            color,
        }
    }

    pub fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        let denom = self.normal.dot(ray.direction);
        
        if denom.abs() < 1e-8 {
            return None; // Ray is parallel to plane
        }

        let t = (self.point - ray.origin).dot(self.normal) / denom;
        
        if t >= t_min && t <= t_max {
            let hit_point = ray.at(t);
            let mut hit_record = HitRecord::new(hit_point, self.normal, t, self.color);
            hit_record.set_face_normal(ray, self.normal);
            Some(hit_record)
        } else {
            None
        }
    }
}