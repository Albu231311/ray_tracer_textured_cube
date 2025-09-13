use crate::vector3::Vector3;
use crate::ray::Ray;

pub struct Camera {
    pub position: Vector3,
    pub target: Vector3,
    pub up: Vector3,
    pub fov: f32,
    pub aspect_ratio: f32,
    pub near: f32,
    pub far: f32,
    
    // Computed values
    forward: Vector3,
    right: Vector3,
    camera_up: Vector3,
    viewport_height: f32,
    viewport_width: f32,
    horizontal: Vector3,
    vertical: Vector3,
    lower_left_corner: Vector3,
}

impl Camera {
    pub fn new(
        position: Vector3,
        target: Vector3,
        up: Vector3,
        fov: f32,
        aspect_ratio: f32,
    ) -> Self {
        let mut camera = Camera {
            position,
            target,
            up,
            fov,
            aspect_ratio,
            near: 0.1,
            far: 100.0,
            forward: Vector3::zero(),
            right: Vector3::zero(),
            camera_up: Vector3::zero(),
            viewport_height: 0.0,
            viewport_width: 0.0,
            horizontal: Vector3::zero(),
            vertical: Vector3::zero(),
            lower_left_corner: Vector3::zero(),
        };
        
        camera.update_vectors();
        camera
    }

    fn update_vectors(&mut self) {
        // Calculate camera coordinate system
        self.forward = (self.target - self.position).normalize();
        self.right = self.forward.cross(self.up).normalize();
        self.camera_up = self.right.cross(self.forward);

        // Calculate viewport dimensions
        let theta = self.fov.to_radians();
        let half_height = (theta / 2.0).tan();
        self.viewport_height = 2.0 * half_height;
        self.viewport_width = self.aspect_ratio * self.viewport_height;

        // Calculate viewport vectors
        self.horizontal = self.right * self.viewport_width;
        self.vertical = self.camera_up * self.viewport_height;
        
        // Calculate lower left corner of viewport - this is key for ray generation
        self.lower_left_corner = self.position 
            + self.forward * 1.0  // Distance to viewport plane
            - self.horizontal * 0.5 
            - self.vertical * 0.5;
    }

    pub fn get_ray(&self, u: f32, v: f32) -> Ray {
        let direction = self.lower_left_corner 
            + self.horizontal * u 
            + self.vertical * v 
            - self.position;
        
        Ray::new(self.position, direction)
    }
}