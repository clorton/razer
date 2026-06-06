use extendr_api::prelude::*;

mod epidemic;
mod distributions;
mod migration;
mod column;
mod bincount;
mod sir;
mod vitals;
mod pyramid;
mod kmestimator;
mod mortality;
mod measles;
mod births;
pub use distributions::Distribution;
pub use column::Column;

extendr_module! {
    mod razer;
    use epidemic;
    use distributions;
    use migration;
    use column;
    use bincount;
    use sir;
    use vitals;
    use pyramid;
    use kmestimator;
    use mortality;
    use measles;
    use births;
}
