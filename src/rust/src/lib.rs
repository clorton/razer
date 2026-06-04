use extendr_api::prelude::*;

mod laser_frame;
mod epidemic;
mod distributions;
pub use laser_frame::LaserFrame;
pub use distributions::Distribution;

extendr_module! {
    mod razer;
    use laser_frame;
    use epidemic;
    use distributions;
}
