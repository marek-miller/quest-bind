#![allow(clippy::cast_sign_loss)]

use super::*;

#[test]
fn create_qureg_01() -> Result<(), QuestError> {
    let env = QuESTEnv::new();
    let _ = Qureg::try_new(1, &env)?;
    let _ = Qureg::try_new(5, &env)?;

    let _ = Qureg::try_new(0, &env).unwrap_err();
    Ok(())
}

#[test]
fn create_density_qureg_01() -> Result<(), QuestError> {
    let env = QuESTEnv::new();
    {
        let _ = Qureg::try_new_density(1, &env)?;
        let _ = Qureg::try_new_density(5, &env)?;

        let _ = Qureg::try_new_density(0, &env).unwrap_err();
    }
    Ok(())
}

#[test]
fn create_clone_qureg_01() -> Result<(), QuestError> {
    let env = QuESTEnv::new();
    {
        let qureg = Qureg::try_new_density(2, &env)?;
        let _ = qureg.clone();
    }
    Ok(())
}

#[test]
fn init_complex_matrix_n_02() -> Result<(), QuestError> {
    let mut m = ComplexMatrixN::try_new(2)?;
    init_complex_matrix_n(
        &mut m,
        &[&[1., 2.], &[3., 4.]],
        &[&[11., 12.], &[13., 14.]],
    )?;

    unsafe {
        let row = &*(*m.0.real).cast::<[&[f64; 2]; 2]>();
        assert_eq!(row, &[&[1., 2.,], &[3., 4.]]);
    }
    unsafe {
        let row = &*(*m.0.imag).cast::<[&[f64; 2]; 2]>();
        assert_eq!(row, &[&[11., 12.], &[13., 14.],]);
    }
    Ok(())
}

#[test]
fn init_complex_matrix_n_03() -> Result<(), QuestError> {
    let mut m = ComplexMatrixN::try_new(3)?;
    init_complex_matrix_n(
        &mut m,
        &[&[1., 2., 3.], &[4., 5., 6.], &[7., 8., 9.]],
        &[&[11., 12., 13.], &[14., 15., 16.], &[17., 18., 19.]],
    )?;

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
    Ok(())
}

#[test]
fn create_pauli_hamil_01() {
    let _ = PauliHamil::try_new(1, 1).unwrap();
    let _ = PauliHamil::try_new(2, 3).unwrap();
    let _ = PauliHamil::try_new(3, 2).unwrap();

    let _ = PauliHamil::try_new(0, 1).unwrap_err();
    let _ = PauliHamil::try_new(-1, 1).unwrap_err();
    let _ = PauliHamil::try_new(1, 0).unwrap_err();
    let _ = PauliHamil::try_new(1, -1).unwrap_err();
    let _ = PauliHamil::try_new(0, 0).unwrap_err();
}

#[test]
fn initialize_pauli_hamil_01() {
    use PauliOpType::*;
    let mut hamil = PauliHamil::try_new(2, 2).unwrap();

    init_pauli_hamil(
        &mut hamil,
        &[0.5, -0.5],
        &[PAULI_X, PAULI_Y, PAULI_I, PAULI_I, PAULI_Z, PAULI_X],
    )
    .unwrap();
}

#[test]
fn get_environment_string_01() {
    let env = QuESTEnv::new();
    let env_str = get_environment_string(&env).unwrap();

    assert!(env_str.contains("CUDA="));
    assert!(env_str.contains("OpenMP="));
    assert!(env_str.contains("MPI="));
    assert!(env_str.contains("threads="));
    assert!(env_str.contains("ranks="));
}

#[test]
fn get_quest_seeds_01() {
    let env = QuESTEnv::new();
    let (seeds, num_seeds) = get_quest_seeds(&env);

    assert!(num_seeds > 0);
    assert_eq!(seeds.len(), num_seeds as usize);
}
