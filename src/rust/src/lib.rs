use extendr_api::prelude::*;

mod epidemic;
mod rng;
mod distributions;
mod migration;
mod column;
mod bincount;
mod transmission;
mod steps;
mod vitals;
mod pyramid;
mod kmestimator;
mod mortality;
mod births;
pub use distributions::Distribution;
pub use column::Column;

extendr_module! {
    mod razer;
    use epidemic;
    use rng;
    use distributions;
    use migration;
    use column;
    use bincount;
    use transmission;
    use steps;
    use vitals;
    use pyramid;
    use kmestimator;
    use mortality;
    use births;
}
