mod vector3;
mod ray;
mod camera;
mod cube;
mod light;
mod raytracer;
mod texture; // Módulo de texturas

use raylib::prelude::*;
use raytracer::RayTracer;

const WINDOW_WIDTH: i32 = 800;
const WINDOW_HEIGHT: i32 = 600;

fn main() {
    let (mut rl, thread) = raylib::init()
        .size(WINDOW_WIDTH, WINDOW_HEIGHT)
        .title("Textured Cube Ray Tracer")
        .build();

    rl.set_target_fps(60);

    // Create raytracer
    let raytracer = RayTracer::new(WINDOW_WIDTH as u32, WINDOW_HEIGHT as u32);
    
    println!("🎮 Starting textured cube ray tracer...");
    let start_time = std::time::Instant::now();
    
    // Render the scene
    let pixel_data = raytracer.render();
    
    let render_time = start_time.elapsed();
    println!("✅ Ray tracing completed in {:.2} seconds", render_time.as_secs_f32());

    // Create image from pixel data
    let image = Image::gen_image_color(WINDOW_WIDTH, WINDOW_HEIGHT, Color::BLACK);
    
    // Copy our rendered pixels to raylib image
    unsafe {
        let image_data = std::slice::from_raw_parts_mut(
            image.data as *mut u8,
            (WINDOW_WIDTH * WINDOW_HEIGHT * 4) as usize,
        );
        image_data.copy_from_slice(&pixel_data);
    }

    // Create texture from image
    let texture = rl.load_texture_from_image(&thread, &image).expect("Failed to create texture");

    // Main game loop
    while !rl.window_should_close() {
        let mut d = rl.begin_drawing(&thread);
        
        d.clear_background(Color::BLACK);
        
        // Draw the ray traced image
        d.draw_texture(&texture, 0, 0, Color::WHITE);
        
        // Draw UI text
        d.draw_text(
            &format!("Textured Cube - Rendered in {:.2}s", render_time.as_secs_f32()),
            10,
            10,
            20,
            Color::WHITE,
        );
        
        d.draw_text(
            "ESC to exit",
            10,
            WINDOW_HEIGHT - 30,
            20,
            Color::WHITE,
        );
    }
    
    println!("🎮 Ray tracer finished!");
}