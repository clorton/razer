use extendr_api::prelude::*;

mod laser_frame;
mod epidemic;
mod distributions;
mod migration;
mod column;
pub use laser_frame::LaserFrame;
pub use distributions::Distribution;
pub use column::Column;

extendr_module! {
    mod razer;
    use laser_frame;
    use epidemic;
    use distributions;
    use migration;
    use column;
}
