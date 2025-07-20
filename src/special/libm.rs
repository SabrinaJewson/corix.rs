unsafe extern "C" {
    /// Evaluates `x * 2 ^ exp`.
    #[link_name = "ldexp"]
    pub safe fn mul_exp2(x: f64, exp: i32) -> f64;
}
