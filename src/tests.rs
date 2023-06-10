#![allow(clippy::cast_sign_loss)]

use super::*;

#[test]
fn init_complex_matrix_n_02() {
    let mut m = create_complex_matrix_n(2);
    init_complex_matrix_n(
        &mut m,
        &[&[1., 2.], &[3., 4.]],
        &[&[11., 12.], &[13., 14.]],
    );

    unsafe {
        let row = &*(*m.0.real).cast::<[&[f64; 2]; 2]>();
        assert_eq!(row, &[&[1., 2.,], &[3., 4.]]);
    }
    unsafe {
        let row = &*(*m.0.imag).cast::<[&[f64; 2]; 2]>();
        assert_eq!(row, &[&[11., 12.], &[13., 14.],]);
    }
    destroy_complex_matrix_n(m);
}

#[test]
fn init_complex_matrix_n_03() {
    let mut m = create_complex_matrix_n(3);
    init_complex_matrix_n(
        &mut m,
        &[&[1., 2., 3.], &[4., 5., 6.], &[7., 8., 9.]],
        &[&[11., 12., 13.], &[14., 15., 16.], &[17., 18., 19.]],
    );

    unsafe {
        let row = &*(*m.0.real).cast::<[&[f64; 3]; 3]>();
        assert_eq!(row, &[&[1., 2., 3.], &[4., 5., 6.], &[7., 8., 9.]]);
    }
    unsafe {
        let row = &*(*m.0.imag).cast::<[&[f64; 3]; 3]>();
        assert_eq!(
            row,
            &[&[11., 12., 13.], &[14., 15., 16.], &[17., 18., 19.]]
        );
    }
    destroy_complex_matrix_n(m);
}

#[test]
fn get_environment_string_01() {
    let env = create_quest_env();
    let env_str = get_environment_string(&env).unwrap();

    assert!(env_str.contains("CUDA="));
    assert!(env_str.contains("OpenMP="));
    assert!(env_str.contains("MPI="));
    assert!(env_str.contains("threads="));
    assert!(env_str.contains("ranks="));

    destroy_quest_env(env);
}

#[test]
fn get_quest_seeds_01() {
    let env = create_quest_env();
    let (seeds, num_seeds) = get_quest_seeds(&env);

    assert!(num_seeds > 0);
    assert_eq!(seeds.len(), num_seeds as usize);

    destroy_quest_env(env);
}
