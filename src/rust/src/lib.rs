use extendr_api::prelude::*;

mod laser_frame;
mod epidemic;
pub use laser_frame::LaserFrame;

extendr_module! {
    mod razer;
    use laser_frame;
    use epidemic;
}
