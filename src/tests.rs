#![allow(clippy::cast_sign_loss)]

use super::*;

#[test]
fn create_qureg_01() -> Result<(), QuestError> {
    let env = &QuestEnv::new();
    let _ = Qureg::try_new(1, env)?;
    let _ = Qureg::try_new(5, env)?;

    let _ = Qureg::try_new(0, env).unwrap_err();
    Ok(())
}

#[test]
fn create_density_qureg_01() -> Result<(), QuestError> {
    let env = &QuestEnv::new();
    {
        let _ = Qureg::try_new_density(1, env)?;
        let _ = Qureg::try_new_density(5, env)?;

        let _ = Qureg::try_new_density(0, env).unwrap_err();
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
fn create_diagonal_op_01() {
    let env = &QuestEnv::new();

    let _ = DiagonalOp::try_new(1, env).unwrap();
    let _ = DiagonalOp::try_new(0, env).unwrap_err();
    let _ = DiagonalOp::try_new(-1, env).unwrap_err();
}

#[test]
fn set_diagonal_op_elems_01() {
    let env = &QuestEnv::new();
    let mut op = DiagonalOp::try_new(3, env).unwrap();

    let num_elems = 3;
    let re = [1., 2., 3.];
    let im = [1., 2., 3.];
    set_diagonal_op_elems(&mut op, 0, &re, &im, num_elems).unwrap();
    set_diagonal_op_elems(&mut op, -1, &re, &im, num_elems).unwrap_err();
    set_diagonal_op_elems(&mut op, 9, &re, &im, 3).unwrap_err();
}

#[test]
fn apply_diagonal_op_01() {
    let env = &QuestEnv::new();
    let mut qureg = Qureg::try_new(2, env).unwrap();
    let mut op = DiagonalOp::try_new(2, env).unwrap();

    init_diagonal_op(&mut op, &[1., 2., 3., 4.], &[5., 6., 7., 8.]).unwrap();
    apply_diagonal_op(&mut qureg, &op).unwrap();

    let mut op = DiagonalOp::try_new(1, env).unwrap();
    init_diagonal_op(&mut op, &[1., 2.], &[5., 6.]).unwrap();
    apply_diagonal_op(&mut qureg, &op).unwrap_err();
}

#[test]
fn calc_expec_diagonal_op_() {
    let env = &QuestEnv::new();
    let mut qureg = Qureg::try_new(2, env).unwrap();
    let mut op = DiagonalOp::try_new(2, env).unwrap();

    init_plus_state(&mut qureg);
    init_diagonal_op(&mut op, &[1., 2., 3., 4.], &[5., 6., 7., 8.]).unwrap();

    let _ = calc_expec_diagonal_op(&qureg, &op).unwrap();

    let mut op = DiagonalOp::try_new(1, env).unwrap();
    init_diagonal_op(&mut op, &[1., 2.], &[5., 6.]).unwrap();
    let _ = calc_expec_diagonal_op(&qureg, &op).unwrap_err();
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
fn set_amps_01() {
    let env = &QuestEnv::new();
    let mut qureg = Qureg::try_new(3, env).unwrap();

    let num_amps = 4;
    let re = [1., 2., 3., 4.];
    let im = [1., 2., 3., 4.];

    set_amps(&mut qureg, 0, &re, &im, num_amps).unwrap();

    assert!((get_real_amp(&qureg, 0).unwrap() - 1.).abs() < f64::EPSILON);

    set_amps(&mut qureg, 9, &re, &im, 4).unwrap_err();
    set_amps(&mut qureg, 7, &re, &im, 4).unwrap_err();
    set_amps(&mut qureg, 3, &re, &im, 9).unwrap_err();
    set_amps(&mut qureg, -1, &re, &im, 9).unwrap_err();
    set_amps(&mut qureg, 1, &re, &im, -9).unwrap_err();
}

#[test]
fn set_density_amps_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new_density(3, env).unwrap();

    let re = &[1., 2., 3., 4.];
    let im = &[1., 2., 3., 4.];

    set_density_amps(qureg, 0, 0, re, im, 4).unwrap();
    assert!(
        (get_density_amp(qureg, 0, 0).unwrap().re - 1.).abs() < f64::EPSILON
    );

    set_amps(qureg, 0, re, im, 4).unwrap_err();

    set_density_amps(qureg, 0, 9, re, im, 4).unwrap_err();
    set_density_amps(qureg, 8, 7, re, im, 4).unwrap_err();
    set_density_amps(qureg, 0, -1, re, im, 9).unwrap_err();
    set_density_amps(qureg, 0, 1, re, im, -9).unwrap_err();
}

#[test]
fn phase_shift_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();

    phase_shift(qureg, 0, 0.0).unwrap();
    phase_shift(qureg, 1, 0.5).unwrap();
    phase_shift(qureg, 2, 1.0).unwrap();

    phase_shift(qureg, 3, 0.0).unwrap_err();
    phase_shift(qureg, -11, 0.0).unwrap_err();
}

#[test]
fn controlled_phase_shift_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();

    controlled_phase_shift(qureg, 0, 1, 0.5).unwrap();
    controlled_phase_shift(qureg, 0, 2, 0.5).unwrap();

    controlled_phase_shift(qureg, 0, 3, 0.5).unwrap_err();
    controlled_phase_shift(qureg, -1, 1, 0.5).unwrap_err();
}

#[test]
fn multi_controlled_phase_shift_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(4, env).unwrap();
    multi_controlled_phase_shift(qureg, &[0, 1, 3], 3, 0.5).unwrap();
    multi_controlled_phase_shift(qureg, &[0, 1, 3], 2, 0.5).unwrap();

    multi_controlled_phase_shift(qureg, &[0, 4, 3, 4], 2, 0.5).unwrap_err();
    multi_controlled_phase_shift(qureg, &[0, 1], 3, 0.5).unwrap_err();
    multi_controlled_phase_shift(qureg, &[0, 1], -1, 0.5).unwrap_err();
}

#[test]
fn controlled_phase_flip_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();

    controlled_phase_flip(qureg, 0, 1).unwrap();
    controlled_phase_flip(qureg, 0, 2).unwrap();

    controlled_phase_flip(qureg, 0, 3).unwrap_err();
    controlled_phase_flip(qureg, -1, 1).unwrap_err();
}

#[test]
fn multi_controlled_phase_flip_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(4, env).unwrap();
    multi_controlled_phase_flip(qureg, &[0, 1, 3]).unwrap();
    multi_controlled_phase_flip(qureg, &[0, 1, 3]).unwrap();

    multi_controlled_phase_flip(qureg, &[0, 4, 3, 4]).unwrap_err();
    multi_controlled_phase_flip(qureg, &[0, 7, -1]).unwrap_err();
}

#[test]
fn s_gate_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    s_gate(qureg, 0).unwrap();
    assert!((get_imag_amp(qureg, 0).unwrap()).abs() < f64::EPSILON);

    pauli_x(qureg, 0).unwrap();
    s_gate(qureg, 0).unwrap();

    let amp = get_imag_amp(qureg, 1).unwrap();
    assert!((amp - 1.).abs() < f64::EPSILON);

    s_gate(qureg, -1).unwrap_err();
    s_gate(qureg, 3).unwrap_err();
}

#[test]
fn t_gate_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    t_gate(qureg, 0).unwrap();
    t_gate(qureg, -1).unwrap_err();
    t_gate(qureg, 3).unwrap_err();
}

#[test]
fn get_amp_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_plus_state(qureg);

    get_amp(qureg, 0).unwrap();
    get_amp(qureg, 1).unwrap();
    get_amp(qureg, 2).unwrap();
    get_amp(qureg, 3).unwrap();

    get_amp(qureg, 4).unwrap_err();
    get_amp(qureg, -1).unwrap_err();
}

#[test]
fn get_real_amp_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_plus_state(qureg);

    get_real_amp(qureg, 0).unwrap();
    get_real_amp(qureg, 1).unwrap();
    get_real_amp(qureg, 2).unwrap();
    get_real_amp(qureg, 3).unwrap();

    get_real_amp(qureg, 4).unwrap_err();
    get_real_amp(qureg, -1).unwrap_err();
}

#[test]
fn get_imag_amp_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_plus_state(qureg);

    get_imag_amp(qureg, 0).unwrap();
    get_imag_amp(qureg, 1).unwrap();
    get_imag_amp(qureg, 2).unwrap();
    get_imag_amp(qureg, 3).unwrap();

    get_imag_amp(qureg, 4).unwrap_err();
    get_imag_amp(qureg, -1).unwrap_err();
}

#[test]
fn get_prob_amp_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_plus_state(qureg);

    get_prob_amp(qureg, 0).unwrap();
    get_prob_amp(qureg, 1).unwrap();
    get_prob_amp(qureg, 2).unwrap();
    get_prob_amp(qureg, 3).unwrap();

    get_prob_amp(qureg, 4).unwrap_err();
    get_prob_amp(qureg, -1).unwrap_err();
}

#[test]
fn get_density_amp_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();

    get_density_amp(qureg, 0, 0).unwrap_err();
    get_density_amp(qureg, 1, 0).unwrap_err();
    get_density_amp(qureg, 2, 0).unwrap_err();
    get_density_amp(qureg, 3, 0).unwrap_err();
    get_density_amp(qureg, -1, 5).unwrap_err();
    get_density_amp(qureg, 4, 0).unwrap_err();
}

#[test]
fn get_density_amp_02() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new_density(2, env).unwrap();

    get_density_amp(qureg, 0, 0).unwrap();
    get_density_amp(qureg, 1, 0).unwrap();
    get_density_amp(qureg, 2, 0).unwrap();
    get_density_amp(qureg, 3, 0).unwrap();
    get_density_amp(qureg, -1, 0).unwrap_err();
    get_density_amp(qureg, 4, 0).unwrap_err();
}

#[test]
fn compact_unitary_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    let norm = std::f64::consts::SQRT_2.recip();
    let alpha = Qcomplex::new(0., norm);
    let beta = Qcomplex::new(0., norm);

    compact_unitary(qureg, 0, alpha, beta).unwrap();
    compact_unitary(qureg, 1, alpha, beta).unwrap();

    compact_unitary(qureg, 4, alpha, beta).unwrap_err();
    compact_unitary(qureg, -1, alpha, beta).unwrap_err();
}

#[test]
fn compact_unitary_02() {
    // env_logger::init();
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    // this doesn't define a unitary matrix
    let alpha = Qcomplex::new(1., 2.);
    let beta = Qcomplex::new(2., 1.);

    compact_unitary(qureg, 0, alpha, beta).unwrap_err();
    compact_unitary(qureg, 1, alpha, beta).unwrap_err();

    compact_unitary(qureg, 4, alpha, beta).unwrap_err();
    compact_unitary(qureg, -1, alpha, beta).unwrap_err();
}

#[test]
fn unitary_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    let norm = std::f64::consts::SQRT_2.recip();
    let mtr = ComplexMatrix2::new(
        [[norm, norm], [norm, -norm]],
        [[0., 0.], [0., 0.]],
    );
    unitary(qureg, 0, &mtr).unwrap();
    unitary(qureg, 1, &mtr).unwrap();
    unitary(qureg, 2, &mtr).unwrap_err();
    unitary(qureg, -1, &mtr).unwrap_err();
}

#[test]
fn unitary_02() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    // This isn't a unitary
    let mtr = ComplexMatrix2::new([[1., 2.], [0., 0.]], [[0., 0.], [1., 2.]]);
    unitary(qureg, 0, &mtr).unwrap_err();
    unitary(qureg, 1, &mtr).unwrap_err();
    unitary(qureg, 2, &mtr).unwrap_err();
    unitary(qureg, -1, &mtr).unwrap_err();
}

#[test]
fn rotate_x_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();
    let theta = std::f64::consts::PI;
    rotate_x(qureg, 0, theta).unwrap();
    rotate_x(qureg, 1, theta).unwrap();
    rotate_x(qureg, 2, theta).unwrap();

    rotate_x(qureg, 3, theta).unwrap_err();
    rotate_x(qureg, -1, theta).unwrap_err();
}

#[test]
fn rotate_y_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();
    let theta = std::f64::consts::PI;
    rotate_y(qureg, 0, theta).unwrap();
    rotate_y(qureg, 1, theta).unwrap();
    rotate_y(qureg, 2, theta).unwrap();

    rotate_y(qureg, 3, theta).unwrap_err();
    rotate_y(qureg, -1, theta).unwrap_err();
}

#[test]
fn rotate_z_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();
    let theta = std::f64::consts::PI;
    rotate_z(qureg, 0, theta).unwrap();
    rotate_z(qureg, 1, theta).unwrap();
    rotate_z(qureg, 2, theta).unwrap();

    rotate_z(qureg, 3, theta).unwrap_err();
    rotate_z(qureg, -1, theta).unwrap_err();
}

#[test]
fn rotate_around_axis_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();
    init_zero_state(qureg);

    let angle = 0.;
    let axis = &Vector::new(0., 0., 1.);
    rotate_around_axis(qureg, 0, angle, axis).unwrap();
    rotate_around_axis(qureg, 1, angle, axis).unwrap();
    rotate_around_axis(qureg, 2, angle, axis).unwrap();

    rotate_around_axis(qureg, 3, angle, axis).unwrap_err();
    rotate_around_axis(qureg, -1, angle, axis).unwrap_err();
}

#[test]
fn rotate_around_axis_02() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();
    init_zero_state(qureg);

    let angle = 0.;
    // zero vector should throw an exception
    let axis = &Vector::new(0., 0., 0.);
    rotate_around_axis(qureg, 0, angle, axis).unwrap_err();
    rotate_around_axis(qureg, 1, angle, axis).unwrap_err();
    rotate_around_axis(qureg, 2, angle, axis).unwrap_err();

    rotate_around_axis(qureg, 3, angle, axis).unwrap_err();
    rotate_around_axis(qureg, -1, angle, axis).unwrap_err();
}

#[test]
fn controlled_rotate_x_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();

    controlled_rotate_x(qureg, 1, 0, 0.5).unwrap();
    controlled_rotate_x(qureg, 1, 2, 0.5).unwrap();

    controlled_rotate_x(qureg, 1, 1, 0.5).unwrap_err();
    controlled_rotate_x(qureg, 2, 2, 0.5).unwrap_err();
    controlled_rotate_x(qureg, -1, 2, 0.5).unwrap_err();
    controlled_rotate_x(qureg, 2, -1, 0.5).unwrap_err();
    controlled_rotate_x(qureg, 0, 4, 0.5).unwrap_err();
    controlled_rotate_x(qureg, 4, 0, 0.5).unwrap_err();
}

#[test]
fn controlled_rotate_y_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();

    controlled_rotate_y(qureg, 1, 0, 0.5).unwrap();
    controlled_rotate_y(qureg, 1, 2, 0.5).unwrap();

    controlled_rotate_y(qureg, 1, 1, 0.5).unwrap_err();
    controlled_rotate_y(qureg, 2, 2, 0.5).unwrap_err();
    controlled_rotate_y(qureg, -1, 2, 0.5).unwrap_err();
    controlled_rotate_y(qureg, 2, -1, 0.5).unwrap_err();
    controlled_rotate_y(qureg, 0, 4, 0.5).unwrap_err();
    controlled_rotate_y(qureg, 4, 0, 0.5).unwrap_err();
}

#[test]
fn controlled_rotate_z_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();

    controlled_rotate_z(qureg, 1, 0, 0.5).unwrap();
    controlled_rotate_z(qureg, 1, 2, 0.5).unwrap();

    controlled_rotate_z(qureg, 1, 1, 0.5).unwrap_err();
    controlled_rotate_z(qureg, 2, 2, 0.5).unwrap_err();
    controlled_rotate_z(qureg, -1, 2, 0.5).unwrap_err();
    controlled_rotate_z(qureg, 2, -1, 0.5).unwrap_err();
    controlled_rotate_z(qureg, 0, 4, 0.5).unwrap_err();
    controlled_rotate_z(qureg, 4, 0, 0.5).unwrap_err();
}

#[test]
fn controlled_rotate_around_axis_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();
    let vector = &Vector::new(0., 0., 1.);

    controlled_rotate_around_axis(qureg, 1, 0, 0.5, vector).unwrap();
    controlled_rotate_around_axis(qureg, 1, 2, 0.5, vector).unwrap();

    controlled_rotate_around_axis(qureg, 1, 1, 0.5, vector).unwrap_err();
    controlled_rotate_around_axis(qureg, 2, 2, 0.5, vector).unwrap_err();
    controlled_rotate_around_axis(qureg, -1, 2, 0.5, vector).unwrap_err();
    controlled_rotate_around_axis(qureg, 2, -1, 0.5, vector).unwrap_err();
    controlled_rotate_around_axis(qureg, 0, 4, 0.5, vector).unwrap_err();
    controlled_rotate_around_axis(qureg, 4, 0, 0.5, vector).unwrap_err();
}

#[test]
fn controlled_rotate_around_axis_02() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();
    // vector cannot be zero
    let vector = &Vector::new(0., 0., 0.);

    controlled_rotate_around_axis(qureg, 1, 0, 0.5, vector).unwrap_err();
    controlled_rotate_around_axis(qureg, 1, 2, 0.5, vector).unwrap_err();
}

#[test]
fn get_quest_seeds_01() {
    let env = &QuestEnv::new();
    let (seeds, num_seeds) = get_quest_seeds(env);

    assert!(num_seeds > 0);
    assert_eq!(seeds.len(), num_seeds as usize);
}
