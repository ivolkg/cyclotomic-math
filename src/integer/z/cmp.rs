//! Implementations to compare [`Z`] with other values.
//! This uses the traits from [`std::cmp`].

use super::Z;

use flint_sys::fmpz::fmpz_equal;

impl PartialEq for Z {
    /// Checks if two integers are equal. Used by the `==` and `!=` operators.
    ///
    /// Parameters:
    /// - other: the other value that is used to compare the elements
    ///
    /// Returns `true` if the elements are equal, otherwise `false`.
    ///
    /// # Example
    /// ```rust
    /// use math::integer::Z;
    /// let a: Z = Z::from(42);
    /// let b: Z = Z::from(24);
    ///
    /// // These are all equivalent and return false.
    /// let compared: bool = (a == b);
    /// # assert!(!compared);
    /// let compared: bool = (&a == &b);
    /// # assert!(!compared);
    /// let compared: bool = (a.eq(&b));
    /// # assert!(!compared);
    /// let compared: bool = (Z::eq(&a,&b));
    /// # assert!(!compared);
    /// ```
    fn eq(&self, other: &Self) -> bool {
        unsafe { 1 == fmpz_equal(&self.value, &other.value) }
    }
}

// With the [`Eq`] trait, `a == a` is always true.
// This is not guaranteed by the [`PartialEq`] trait.
impl Eq for Z {}

#[cfg(test)]
mod test_equal {
    use super::Z;

    // Assert 42 != 24
    #[test]
    fn not_equal() {
        let a = Z::from_i64(42);
        let b = Z::from_i64(24);

        assert!(a != b);
        assert!(&a != &b);
        assert_ne!(a, b);
        assert_ne!(b, a);
    }

    // Test for equality with a large number that uses FLINT's pointer representation.
    #[test]
    fn equal_u64_max() {
        let a = Z::from_u64(u64::MAX);
        let b = Z::from_u64(u64::MAX);

        assert!(a == b);
        assert!(&a == &b);
        assert_eq!(a, b);
        assert_eq!(b, a);
        assert_eq!(a, a);
    }

    // Test for non equality with large numbers that use FLINT's pointer representation.
    #[test]
    fn not_equal_large() {
        let a = Z::from_u64(u64::MAX);
        let b = Z::from_u64(u64::MAX - 1);

        assert_ne!(a, b);
        assert_ne!(&a, &b);
        assert_ne!(b, a);
    }
}
