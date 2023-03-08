//! This module contains basic traits for this library. These include
//! specific traits for matrices and polynomials.

/// Is implemented by polynomials to evaluate it for a certain input.
pub trait Evaluate<U, V> {
    /// Evaluates the object, e.g. polynomial or a matrix of polynomials,
    /// for a given input value.
    ///
    /// Parameters:
    /// - `value`: The value with which to evaluate the object.
    ///
    /// Returns the evaluation of the object.
    fn evaluate<T: Into<U>>(&self, value: T) -> V;
}