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
fn get_matrix_n_elem_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();

    init_zero_state(qureg);
    let amp = get_imag_amp(qureg, 0).unwrap();
    assert_eq!(amp, 0.);

    let mtr = &mut ComplexMatrixN::try_new(1).unwrap();
    init_complex_matrix_n(
        mtr,
        &[&[1., 2.], &[3., 4.]],
        &[&[11., 12.], &[13., 14.]],
    )
    .unwrap();
    apply_matrix_n(qureg, &[0], 1, mtr).unwrap();
    let amp = get_imag_amp(qureg, 0).unwrap();
    assert!((amp - 11.).abs() < EPSILON);
}

#[test]
fn init_complex_matrix_n_02() {
    let mut m = ComplexMatrixN::try_new(1).unwrap();
    init_complex_matrix_n(
        &mut m,
        &[&[1., 2.], &[3., 4.]],
        &[&[11., 12.], &[13., 14.]],
    )
    .unwrap();

    unsafe {
        let elem_ptr = (*(m.0.real.add(0))).add(0);
        assert_eq!(*elem_ptr, 1.);
        let elem_ptr = (*(m.0.real.add(0))).add(1);
        assert_eq!(*elem_ptr, 2.);
        let elem_ptr = (*(m.0.real.add(1))).add(0);
        assert_eq!(*elem_ptr, 3.);
        let elem_ptr = (*(m.0.real.add(1))).add(1);
        assert_eq!(*elem_ptr, 4.);
    }

    unsafe {
        let elem_ptr = (*(m.0.imag.add(0))).add(0);
        assert_eq!(*elem_ptr, 11.);
        let elem_ptr = (*(m.0.imag.add(0))).add(1);
        assert_eq!(*elem_ptr, 12.);
        let elem_ptr = (*(m.0.imag.add(1))).add(0);
        assert_eq!(*elem_ptr, 13.);
        let elem_ptr = (*(m.0.imag.add(1))).add(1);
        assert_eq!(*elem_ptr, 14.);
    }
}

#[test]
fn init_complex_matrix_n_03() {
    let mut m = ComplexMatrixN::try_new(2).unwrap();
    init_complex_matrix_n(
        &mut m,
        &[
            &[111., 112., 113., 114.],
            &[115., 116., 117., 118.],
            &[119., 120., 121., 122.],
            &[123., 124., 125., 126.],
        ],
        &[
            &[211., 212., 213., 214.],
            &[215., 216., 217., 218.],
            &[219., 220., 221., 222.],
            &[223., 224., 225., 226.],
        ],
    )
    .unwrap();

    unsafe {
        let elem_ptr = (*(m.0.real.add(0))).add(0);
        assert_eq!(*elem_ptr, 111.);
        let elem_ptr = (*(m.0.real.add(0))).add(1);
        assert_eq!(*elem_ptr, 112.);
        let elem_ptr = (*(m.0.real.add(0))).add(2);
        assert_eq!(*elem_ptr, 113.);
        let elem_ptr = (*(m.0.real.add(0))).add(3);
        assert_eq!(*elem_ptr, 114.);
    }
    unsafe {
        let elem_ptr = (*(m.0.real.add(1))).add(0);
        assert_eq!(*elem_ptr, 115.);
        let elem_ptr = (*(m.0.real.add(1))).add(1);
        assert_eq!(*elem_ptr, 116.);
        let elem_ptr = (*(m.0.real.add(1))).add(2);
        assert_eq!(*elem_ptr, 117.);
        let elem_ptr = (*(m.0.real.add(1))).add(3);
        assert_eq!(*elem_ptr, 118.);
    }

    unsafe {
        let elem_ptr = (*(m.0.imag.add(0))).add(0);
        assert_eq!(*elem_ptr, 211.);
        let elem_ptr = (*(m.0.imag.add(0))).add(1);
        assert_eq!(*elem_ptr, 212.);
        let elem_ptr = (*(m.0.imag.add(0))).add(2);
        assert_eq!(*elem_ptr, 213.);
        let elem_ptr = (*(m.0.imag.add(0))).add(3);
        assert_eq!(*elem_ptr, 214.);
    }
    unsafe {
        let elem_ptr = (*(m.0.imag.add(1))).add(0);
        assert_eq!(*elem_ptr, 215.);
        let elem_ptr = (*(m.0.imag.add(1))).add(1);
        assert_eq!(*elem_ptr, 216.);
        let elem_ptr = (*(m.0.imag.add(1))).add(2);
        assert_eq!(*elem_ptr, 217.);
        let elem_ptr = (*(m.0.imag.add(1))).add(3);
        assert_eq!(*elem_ptr, 218.);
    }
}

#[test]
fn complex_matrix_n_row_slice_02() {
    let mtr = &mut ComplexMatrixN::try_new(2).unwrap();

    let row = mtr.row_real_as_mut_slice(0);
    row[0] = 1.;
    row[1] = 2.;
    assert_eq!(row.len(), 4);

    let row = mtr.row_real_as_slice(0);
    assert_eq!(row.len(), 4);
    assert_eq!(row[0], 1.);
    assert_eq!(row[1], 2.);

    let row = mtr.row_imag_as_mut_slice(0);
    row[0] = 3.;
    row[1] = 4.;
    assert_eq!(row.len(), 4);

    let row = mtr.row_imag_as_slice(0);
    assert_eq!(row.len(), 4);
    assert_eq!(row[0], 3.);
    assert_eq!(row[1], 4.);
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

    assert!((get_real_amp(&qureg, 0).unwrap() - 1.).abs() < EPSILON);

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
    assert!((get_density_amp(qureg, 0, 0).unwrap().re - 1.).abs() < EPSILON);

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
    assert!((get_imag_amp(qureg, 0).unwrap()).abs() < EPSILON);

    pauli_x(qureg, 0).unwrap();
    s_gate(qureg, 0).unwrap();

    let amp = get_imag_amp(qureg, 1).unwrap();
    assert!((amp - 1.).abs() < EPSILON);

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

    let norm = SQRT_2.recip();
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

    let norm = SQRT_2.recip();
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
    let theta = PI;
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
    let theta = PI;
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
    let theta = PI;
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
fn controlled_compact_unitary_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    let norm = SQRT_2.recip();
    let alpha = Qcomplex::new(0., norm);
    let beta = Qcomplex::new(0., norm);

    controlled_compact_unitary(qureg, 0, 1, alpha, beta).unwrap();
    controlled_compact_unitary(qureg, 1, 0, alpha, beta).unwrap();

    controlled_compact_unitary(qureg, 1, 1, alpha, beta).unwrap_err();
    controlled_compact_unitary(qureg, 2, 2, alpha, beta).unwrap_err();
    controlled_compact_unitary(qureg, 4, 1, alpha, beta).unwrap_err();
    controlled_compact_unitary(qureg, -1, 1, alpha, beta).unwrap_err();
}

#[test]
fn controlled_compact_unitary_02() {
    // env_logger::init();
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    // this doesn't define a unitary matrix
    let alpha = Qcomplex::new(1., 2.);
    let beta = Qcomplex::new(2., 1.);

    controlled_compact_unitary(qureg, 0, 1, alpha, beta).unwrap_err();
    controlled_compact_unitary(qureg, 1, 0, alpha, beta).unwrap_err();

    controlled_compact_unitary(qureg, 1, 1, alpha, beta).unwrap_err();
    controlled_compact_unitary(qureg, 4, 1, alpha, beta).unwrap_err();
    controlled_compact_unitary(qureg, -1, 2, alpha, beta).unwrap_err();
}

#[test]
fn controlled_unitary_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    let norm = SQRT_2.recip();
    let mtr = &ComplexMatrix2::new(
        [[norm, norm], [norm, -norm]],
        [[0., 0.], [0., 0.]],
    );

    controlled_unitary(qureg, 0, 1, mtr).unwrap();
    controlled_unitary(qureg, 1, 0, mtr).unwrap();

    controlled_unitary(qureg, 1, 1, mtr).unwrap_err();
    controlled_unitary(qureg, 2, 2, mtr).unwrap_err();
    controlled_unitary(qureg, 4, 1, mtr).unwrap_err();
    controlled_unitary(qureg, -1, 1, mtr).unwrap_err();
}

#[test]
fn controlled_unitary_02() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    // this doesn't define a unitary matrix
    let mtr = &ComplexMatrix2::new([[1., 2.], [3., 4.]], [[5., 6.], [7., 8.]]);

    controlled_unitary(qureg, 0, 1, mtr).unwrap_err();
    controlled_unitary(qureg, 1, 0, mtr).unwrap_err();

    controlled_unitary(qureg, 1, 1, mtr).unwrap_err();
    controlled_unitary(qureg, 2, 2, mtr).unwrap_err();
    controlled_unitary(qureg, 4, 1, mtr).unwrap_err();
    controlled_unitary(qureg, -1, 1, mtr).unwrap_err();
}

#[test]
fn multi_controlled_unitary_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();
    init_zero_state(qureg);

    let norm = SQRT_2.recip();
    let mtr = &ComplexMatrix2::new(
        [[norm, norm], [norm, -norm]],
        [[0., 0.], [0., 0.]],
    );

    multi_controlled_unitary(qureg, &[0], 2, mtr).unwrap();
    multi_controlled_unitary(qureg, &[0, 1], 2, mtr).unwrap();
    multi_controlled_unitary(qureg, &[1, 0], 2, mtr).unwrap();
    multi_controlled_unitary(qureg, &[1, 2], 0, mtr).unwrap();

    multi_controlled_unitary(qureg, &[1, 1], 1, mtr).unwrap_err();
    multi_controlled_unitary(qureg, &[1, 1], 4, mtr).unwrap_err();
    multi_controlled_unitary(qureg, &[-1, 1], 0, mtr).unwrap_err();
    multi_controlled_unitary(qureg, &[1, 1, 1], 0, mtr).unwrap_err();
    multi_controlled_unitary(qureg, &[0, 1, 2, 3], 0, mtr).unwrap_err();
}

#[test]
fn multi_controlled_unitary_02() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();
    init_zero_state(qureg);

    // this doesn't define a unitary matrix
    let mtr = &ComplexMatrix2::new([[1., 2.], [3., 4.]], [[5., 6.], [7., 8.]]);

    multi_controlled_unitary(qureg, &[0, 1], 2, mtr).unwrap_err();
    multi_controlled_unitary(qureg, &[1, 2], 0, mtr).unwrap_err();
    multi_controlled_unitary(qureg, &[1, 1], 1, mtr).unwrap_err();
    multi_controlled_unitary(qureg, &[1, 1], 4, mtr).unwrap_err();
}

#[test]
fn pauli_x_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    pauli_x(qureg, 0).unwrap();
    pauli_x(qureg, 1).unwrap();

    pauli_x(qureg, 2).unwrap_err();
    pauli_x(qureg, -1).unwrap_err();
}

#[test]
fn pauli_y_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    pauli_y(qureg, 0).unwrap();
    pauli_y(qureg, 1).unwrap();

    pauli_y(qureg, 2).unwrap_err();
    pauli_y(qureg, -1).unwrap_err();
}

#[test]
fn pauli_z_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    pauli_z(qureg, 0).unwrap();
    pauli_z(qureg, 1).unwrap();

    pauli_z(qureg, 2).unwrap_err();
    pauli_z(qureg, -1).unwrap_err();
}

#[test]
fn hadamard_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    hadamard(qureg, 0).unwrap();
    hadamard(qureg, 1).unwrap();
    hadamard(qureg, 2).unwrap_err();
    hadamard(qureg, -1).unwrap_err();
}

#[test]
fn controlled_not_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);
    pauli_x(qureg, 1).unwrap();

    controlled_not(qureg, 1, 0).unwrap();
    controlled_not(qureg, 0, 1).unwrap();

    controlled_not(qureg, 0, 0).unwrap_err();
    controlled_not(qureg, 1, 1).unwrap_err();
    controlled_not(qureg, 1, 2).unwrap_err();
    controlled_not(qureg, 2, 4).unwrap_err();
    controlled_not(qureg, 2, -1).unwrap_err();
}

#[test]
fn multi_qubit_not_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    multi_qubit_not(qureg, &[0, 1]).unwrap();
    multi_qubit_not(qureg, &[1, 0]).unwrap();
    multi_qubit_not(qureg, &[0, 0]).unwrap_err();
    multi_qubit_not(qureg, &[1, 1]).unwrap_err();
    multi_qubit_not(qureg, &[4, 1]).unwrap_err();
    multi_qubit_not(qureg, &[0, -1]).unwrap_err();
}

#[test]
fn multi_controlled_multi_qubit_not_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(4, env).unwrap();
    init_zero_state(qureg);

    multi_controlled_multi_qubit_not(qureg, &[0, 1], &[2, 3]).unwrap();
    multi_controlled_multi_qubit_not(qureg, &[1, 0], &[3, 2]).unwrap();
    multi_controlled_multi_qubit_not(qureg, &[1, 0], &[3]).unwrap();
    multi_controlled_multi_qubit_not(qureg, &[1], &[3, 0]).unwrap();

    multi_controlled_multi_qubit_not(qureg, &[1, 0], &[0]).unwrap_err();
    multi_controlled_multi_qubit_not(qureg, &[0, 0], &[1]).unwrap_err();
    multi_controlled_multi_qubit_not(qureg, &[0, 0], &[-1]).unwrap_err();
    multi_controlled_multi_qubit_not(qureg, &[4, 1], &[0]).unwrap_err();
    multi_controlled_multi_qubit_not(qureg, &[0, 1], &[4]).unwrap_err();
}

#[test]
fn controlled_pauli_y_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    controlled_pauli_y(qureg, 1, 0).unwrap();
    controlled_pauli_y(qureg, 0, 1).unwrap();

    controlled_pauli_y(qureg, 0, 0).unwrap_err();
    controlled_pauli_y(qureg, 1, 1).unwrap_err();
    controlled_pauli_y(qureg, 1, 2).unwrap_err();
    controlled_pauli_y(qureg, 2, 4).unwrap_err();
    controlled_pauli_y(qureg, 2, -1).unwrap_err();
}

#[test]
fn calc_prob_of_outcome_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    let _ = calc_prob_of_outcome(qureg, 0, 0).unwrap();
    let _ = calc_prob_of_outcome(qureg, 0, 1).unwrap();
    let _ = calc_prob_of_outcome(qureg, 1, 0).unwrap();
    let _ = calc_prob_of_outcome(qureg, 1, 1).unwrap();

    let _ = calc_prob_of_outcome(qureg, 0, 2).unwrap_err();
    let _ = calc_prob_of_outcome(qureg, 0, -2).unwrap_err();
    let _ = calc_prob_of_outcome(qureg, 1, 3).unwrap_err();
    let _ = calc_prob_of_outcome(qureg, 4, 0).unwrap_err();
}

#[test]
fn calc_prob_of_all_outcomes_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();
    init_zero_state(qureg);

    let outcome_probs = &mut vec![0.; 4];
    calc_prob_of_all_outcomes(outcome_probs, qureg, &[1, 2]).unwrap();
    calc_prob_of_all_outcomes(outcome_probs, qureg, &[0, 1]).unwrap();
    calc_prob_of_all_outcomes(outcome_probs, qureg, &[0, 2]).unwrap();

    calc_prob_of_all_outcomes(outcome_probs, qureg, &[1, 2, 3]).unwrap_err();
    calc_prob_of_all_outcomes(outcome_probs, qureg, &[0, 0]).unwrap_err();
    calc_prob_of_all_outcomes(outcome_probs, qureg, &[4, 0]).unwrap_err();
    calc_prob_of_all_outcomes(outcome_probs, qureg, &[0, -1]).unwrap_err();
}

#[test]
fn collapse_to_outcome_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();

    init_zero_state(qureg);
    collapse_to_outcome(qureg, 0, 0).unwrap();

    init_zero_state(qureg);
    collapse_to_outcome(qureg, 0, 1).unwrap_err();

    init_zero_state(qureg);
    collapse_to_outcome(qureg, -1, 0).unwrap_err();
    collapse_to_outcome(qureg, 3, 0).unwrap_err();
    collapse_to_outcome(qureg, 1, 3).unwrap_err();
    collapse_to_outcome(qureg, 4, 3).unwrap_err();
}

#[test]
fn measure_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();

    init_zero_state(qureg);

    let _ = measure(qureg, 0).unwrap();
    let _ = measure(qureg, 1).unwrap();
    let _ = measure(qureg, -1).unwrap_err();
    let _ = measure(qureg, 3).unwrap_err();
}

#[test]
fn measure_with_stats_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();

    // Prepare a triplet state `|00> + |11>`
    init_zero_state(qureg);

    let prob = &mut -1.;
    let _ = measure_with_stats(qureg, 0, prob).unwrap();
    let _ = measure_with_stats(qureg, 1, prob).unwrap();
    let _ = measure_with_stats(qureg, -1, prob).unwrap_err();
    let _ = measure_with_stats(qureg, 3, prob).unwrap_err();
}

#[test]
fn calc_inner_product_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);
    let other_qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(other_qureg);

    let _ = calc_inner_product(qureg, other_qureg).unwrap();
}

#[test]
fn calc_inner_product_02() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);
    let other_qureg = &mut Qureg::try_new(1, env).unwrap();
    init_zero_state(other_qureg);

    let _ = calc_inner_product(qureg, other_qureg).unwrap_err();
}

#[test]
fn calc_inner_product_03() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);
    let other_qureg = &mut Qureg::try_new_density(2, env).unwrap();
    init_zero_state(other_qureg);

    let _ = calc_inner_product(qureg, other_qureg).unwrap_err();
}

#[test]
fn calc_density_inner_product_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new_density(2, env).unwrap();
    init_zero_state(qureg);
    let other_qureg = &mut Qureg::try_new_density(2, env).unwrap();
    init_zero_state(other_qureg);

    let _ = calc_density_inner_product(qureg, other_qureg).unwrap();
}

#[test]
fn calc_density_inner_product_02() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new_density(2, env).unwrap();
    init_zero_state(qureg);
    let other_qureg = &mut Qureg::try_new_density(1, env).unwrap();
    init_zero_state(other_qureg);

    let _ = calc_density_inner_product(qureg, other_qureg).unwrap_err();
}

#[test]
fn calc_density_inner_product_03() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new_density(2, env).unwrap();
    init_zero_state(qureg);
    let other_qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(other_qureg);

    let _ = calc_density_inner_product(qureg, other_qureg).unwrap_err();
}

#[test]
fn get_quest_seeds_01() {
    let env = &QuestEnv::new();
    let seeds = get_quest_seeds(env);

    assert!(!seeds.is_empty());
}

#[test]
fn get_quest_seeds_02() {
    let env = &mut QuestEnv::new();
    let seed_array = &[0, 1, 2, 3];
    seed_quest(env, seed_array);
    let seeds = get_quest_seeds(env);

    assert!(!seeds.is_empty());
    assert_eq!(seed_array, seeds);
}

#[test]
fn start_recording_qasm_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();

    start_recording_qasm(qureg);
    hadamard(qureg, 0).and(controlled_not(qureg, 0, 1)).unwrap();
    stop_recording_qasm(qureg);

    print_recorded_qasm(qureg);
}

#[test]
fn mix_dephasing_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new_density(2, env).unwrap();
    init_plus_state(qureg);

    mix_dephasing(qureg, 0, 0.5).unwrap();
    mix_dephasing(qureg, 1, 0.0).unwrap();

    mix_dephasing(qureg, 0, 0.75).unwrap_err();
    mix_dephasing(qureg, 2, 0.25).unwrap_err();
    mix_dephasing(qureg, 0, -0.25).unwrap_err();
    mix_dephasing(qureg, -10, 0.25).unwrap_err();
}

#[test]
fn mix_dephasing_02() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_plus_state(qureg);

    // qureg is not a density matrix
    mix_dephasing(qureg, 0, 0.5).unwrap_err();
    mix_dephasing(qureg, 1, 0.0).unwrap_err();
}

#[test]
fn mix_two_qubit_dephasing_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new_density(3, env).unwrap();
    init_plus_state(qureg);

    mix_two_qubit_dephasing(qureg, 0, 1, 0.75).unwrap();
    mix_two_qubit_dephasing(qureg, 0, 2, 0.75).unwrap();
    mix_two_qubit_dephasing(qureg, 1, 2, 0.75).unwrap();
    mix_two_qubit_dephasing(qureg, 1, 0, 0.75).unwrap();
    mix_two_qubit_dephasing(qureg, 2, 1, 0.75).unwrap();

    mix_two_qubit_dephasing(qureg, 0, 1, 0.99).unwrap_err();
    mix_two_qubit_dephasing(qureg, 2, 1, 0.99).unwrap_err();

    mix_two_qubit_dephasing(qureg, 4, 0, 0.1).unwrap_err();
    mix_two_qubit_dephasing(qureg, 0, 4, 0.1).unwrap_err();

    mix_two_qubit_dephasing(qureg, -1, 0, 0.1).unwrap_err();
    mix_two_qubit_dephasing(qureg, 0, -1, 0.1).unwrap_err();
}

#[test]
fn mix_two_qubit_dephasing_02() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();
    init_plus_state(qureg);

    // qureg is not a density matrix
    mix_two_qubit_dephasing(qureg, 0, 1, 0.75).unwrap_err();
    mix_two_qubit_dephasing(qureg, 0, 2, 0.75).unwrap_err();
}

#[test]
fn mix_depolarising_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new_density(2, env).unwrap();
    init_zero_state(qureg);

    mix_depolarising(qureg, 0, 0.00).unwrap();
    mix_depolarising(qureg, 0, 0.75).unwrap();
    mix_depolarising(qureg, 1, 0.75).unwrap();

    mix_depolarising(qureg, 0, 0.99).unwrap_err();
    mix_depolarising(qureg, 1, 0.99).unwrap_err();
    mix_depolarising(qureg, 0, -0.99).unwrap_err();
    mix_depolarising(qureg, -1, 0.99).unwrap_err();
    mix_depolarising(qureg, -1, -0.99).unwrap_err();
}

#[test]
fn mix_damping_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new_density(2, env).unwrap();
    init_plus_state(qureg);

    mix_damping(qureg, 0, 1.).unwrap();
    mix_damping(qureg, 0, 0.).unwrap();
    mix_damping(qureg, 1, 1.).unwrap();
    mix_damping(qureg, 1, 0.).unwrap();

    mix_damping(qureg, 0, 10.).unwrap_err();
    mix_damping(qureg, 0, -10.).unwrap_err();
    mix_damping(qureg, 1, 10.).unwrap_err();
    mix_damping(qureg, 1, -10.).unwrap_err();
    mix_damping(qureg, 3, 0.5).unwrap_err();
    mix_damping(qureg, -3, 0.5).unwrap_err();
}

#[test]
fn mix_damping_02() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_plus_state(qureg);

    // qureg is not a density matrix
    mix_damping(qureg, 1, 0.).unwrap_err();
    mix_damping(qureg, 0, 0.).unwrap_err();
    // QuEST seg faults here:
    mix_damping(qureg, 0, 0.5).unwrap_err();
}

#[test]
fn mix_two_qubit_depolarising_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new_density(3, env).unwrap();
    init_plus_state(qureg);

    mix_two_qubit_depolarising(qureg, 0, 1, 15. / 16.).unwrap();
    mix_two_qubit_depolarising(qureg, 0, 2, 15. / 16.).unwrap();
    mix_two_qubit_depolarising(qureg, 1, 2, 15. / 16.).unwrap();
    mix_two_qubit_depolarising(qureg, 1, 0, 15. / 16.).unwrap();
    mix_two_qubit_depolarising(qureg, 2, 1, 15. / 16.).unwrap();

    mix_two_qubit_depolarising(qureg, 0, 1, 0.99).unwrap_err();
    mix_two_qubit_depolarising(qureg, 2, 1, 0.99).unwrap_err();

    mix_two_qubit_depolarising(qureg, 4, 0, 0.1).unwrap_err();
    mix_two_qubit_depolarising(qureg, 0, 4, 0.1).unwrap_err();

    mix_two_qubit_depolarising(qureg, -1, 0, 0.1).unwrap_err();
    mix_two_qubit_depolarising(qureg, 0, -1, 0.1).unwrap_err();
}

#[test]
fn mix_two_qubit_depolarising_02() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();
    init_plus_state(qureg);

    // qureg is not a density matrix
    mix_two_qubit_depolarising(qureg, 0, 1, 0.75).unwrap_err();
    mix_two_qubit_depolarising(qureg, 0, 2, 0.75).unwrap_err();
}

#[test]
fn mix_pauli_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new_density(2, env).unwrap();
    init_zero_state(qureg);

    let (prob_x, prob_y, prob_z) = (0.25, 0.25, 0.25);
    mix_pauli(qureg, 0, prob_x, prob_y, prob_z).unwrap();
    mix_pauli(qureg, 1, prob_x, prob_y, prob_z).unwrap();

    mix_pauli(qureg, 2, prob_x, prob_y, prob_z).unwrap_err();
    mix_pauli(qureg, -2, prob_x, prob_y, prob_z).unwrap_err();
}

#[test]
fn mix_pauli_02() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new_density(2, env).unwrap();
    init_zero_state(qureg);

    // this is not a prob distribution
    let (prob_x, prob_y, prob_z) = (0.5, 0.5, 0.5);
    mix_pauli(qureg, 0, prob_x, prob_y, prob_z).unwrap_err();
    mix_pauli(qureg, 1, prob_x, prob_y, prob_z).unwrap_err();
}

#[test]
fn mix_pauli_03() {
    let env = &QuestEnv::new();
    // not a density matrix
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    let (prob_x, prob_y, prob_z) = (0.25, 0.25, 0.25);
    mix_pauli(qureg, 0, prob_x, prob_y, prob_z).unwrap_err();
    mix_pauli(qureg, 1, prob_x, prob_y, prob_z).unwrap_err();
}

#[test]
fn mix_density_matrix_01() {
    let env = &QuestEnv::new();
    let combine_qureg = &mut Qureg::try_new_density(2, env).unwrap();
    let other_qureg = &mut Qureg::try_new_density(2, env).unwrap();

    init_zero_state(combine_qureg);
    init_zero_state(other_qureg);

    mix_density_matrix(combine_qureg, 0.0, other_qureg).unwrap();
    mix_density_matrix(combine_qureg, 0.5, other_qureg).unwrap();
    mix_density_matrix(combine_qureg, 0.99, other_qureg).unwrap();

    mix_density_matrix(combine_qureg, 1.01, other_qureg).unwrap_err();
    mix_density_matrix(combine_qureg, -1.01, other_qureg).unwrap_err();
}

#[test]
fn mix_density_matrix_02() {
    let env = &QuestEnv::new();
    // this is not a density matrix
    let combine_qureg = &mut Qureg::try_new(2, env).unwrap();
    let other_qureg = &mut Qureg::try_new_density(2, env).unwrap();

    init_zero_state(combine_qureg);
    init_zero_state(other_qureg);

    mix_density_matrix(combine_qureg, 0.0, other_qureg).unwrap_err();
}

#[test]
fn mix_density_matrix_03() {
    let env = &QuestEnv::new();
    let combine_qureg = &mut Qureg::try_new_density(2, env).unwrap();
    // this is not a density matrix
    let other_qureg = &mut Qureg::try_new(2, env).unwrap();

    init_zero_state(combine_qureg);
    init_zero_state(other_qureg);

    mix_density_matrix(combine_qureg, 0.0, other_qureg).unwrap_err();
}

#[test]
fn mix_density_matrix_04() {
    let env = &QuestEnv::new();
    // dimensions don't match
    let combine_qureg = &mut Qureg::try_new_density(2, env).unwrap();
    let other_qureg = &mut Qureg::try_new_density(3, env).unwrap();

    init_zero_state(combine_qureg);
    init_zero_state(other_qureg);

    mix_density_matrix(combine_qureg, 0.0, other_qureg).unwrap_err();
}

#[test]
fn calc_purity_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new_density(2, env).unwrap();
    init_zero_state(qureg);

    let _ = calc_purity(qureg).unwrap();

    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);
    let _ = calc_purity(qureg).unwrap_err();
}

#[test]
fn calc_fidelity_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new_density(2, env).unwrap();
    let pure_state = &mut Qureg::try_new(2, env).unwrap();

    init_zero_state(qureg);
    init_zero_state(pure_state);

    let _ = calc_fidelity(qureg, pure_state).unwrap();
}

#[test]
fn calc_fidelity_02() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new_density(3, env).unwrap();
    let pure_state = &mut Qureg::try_new(2, env).unwrap();

    init_zero_state(qureg);
    init_zero_state(pure_state);

    let _ = calc_fidelity(qureg, pure_state).unwrap_err();
}

#[test]
fn calc_fidelity_03() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    let pure_state = &mut Qureg::try_new_density(2, env).unwrap();

    init_zero_state(qureg);
    init_zero_state(pure_state);

    let _ = calc_fidelity(qureg, pure_state).unwrap_err();
}

#[test]
fn swap_gate_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    swap_gate(qureg, 0, 1).unwrap();
    swap_gate(qureg, 1, 0).unwrap();

    swap_gate(qureg, 0, 0).unwrap_err();
    swap_gate(qureg, 1, 1).unwrap_err();
}

#[test]
fn swap_gate_02() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    swap_gate(qureg, 0, 2).unwrap_err();
    swap_gate(qureg, 2, 0).unwrap_err();

    swap_gate(qureg, -1, 0).unwrap_err();
    swap_gate(qureg, 0, -1).unwrap_err();

    swap_gate(qureg, 4, 0).unwrap_err();
    swap_gate(qureg, 0, 4).unwrap_err();
    swap_gate(qureg, 4, 4).unwrap_err();
}

#[test]
fn swap_gate_03() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    // QuEST seg faults here
    swap_gate(qureg, 4, -4).unwrap_err();
    swap_gate(qureg, -4, 4).unwrap_err();
    swap_gate(qureg, -4, -4).unwrap_err();
}

#[test]
fn sqrt_swap_gate_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    sqrt_swap_gate(qureg, 0, 1).unwrap();
    sqrt_swap_gate(qureg, 1, 0).unwrap();

    sqrt_swap_gate(qureg, 0, 0).unwrap_err();
    sqrt_swap_gate(qureg, 1, 1).unwrap_err();
}

#[test]
fn sqrt_swap_gate_02() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    sqrt_swap_gate(qureg, 0, 2).unwrap_err();
    sqrt_swap_gate(qureg, 2, 0).unwrap_err();

    sqrt_swap_gate(qureg, -1, 0).unwrap_err();
    sqrt_swap_gate(qureg, 0, -1).unwrap_err();

    sqrt_swap_gate(qureg, 4, 0).unwrap_err();
    sqrt_swap_gate(qureg, 0, 4).unwrap_err();
    sqrt_swap_gate(qureg, 4, 4).unwrap_err();
}

#[test]
fn sqrt_swap_gate_03() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    // QuEST seg faults here
    sqrt_swap_gate(qureg, 4, -4).unwrap_err();
    sqrt_swap_gate(qureg, -4, 4).unwrap_err();
    sqrt_swap_gate(qureg, -4, -4).unwrap_err();
}

#[test]
fn multi_rotate_z_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    multi_rotate_z(qureg, &[0, 1], PI).unwrap();
    multi_rotate_z(qureg, &[0, 1, 2], PI).unwrap_err();
    multi_rotate_z(qureg, &[0, 2], PI).unwrap_err();
    multi_rotate_z(qureg, &[0, 0], PI).unwrap_err();
}

#[test]
fn multi_state_controlled_unitary_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();
    init_zero_state(qureg);

    let u = &ComplexMatrix2::new([[0., 1.], [1., 0.]], [[0., 0.], [0., 0.]]);
    multi_state_controlled_unitary(qureg, &[1, 2], &[0, 0], 0, u).unwrap();
    multi_state_controlled_unitary(qureg, &[0, 2], &[0, 0], 1, u).unwrap();
    multi_state_controlled_unitary(qureg, &[0, 1], &[0, 0], 2, u).unwrap();

    multi_state_controlled_unitary(qureg, &[0, 1], &[0, 0], 0, u).unwrap_err();
    multi_state_controlled_unitary(qureg, &[0, 1], &[0, 0], 1, u).unwrap_err();
    multi_state_controlled_unitary(qureg, &[0, 0], &[0, 0], 1, u).unwrap_err();

    multi_state_controlled_unitary(qureg, &[0, 0], &[0, 0], 3, u).unwrap_err();
    multi_state_controlled_unitary(qureg, &[0, 0], &[0, 0], -1, u).unwrap_err();
    multi_state_controlled_unitary(qureg, &[4, 0], &[0, 0], 1, u).unwrap_err();
    multi_state_controlled_unitary(qureg, &[4, -1], &[0, 0], 1, u).unwrap_err();
}

#[test]
fn apply_matrix2_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    let u = &ComplexMatrix2::new([[0., 1.], [1., 0.]], [[0., 0.], [0., 0.]]);

    apply_matrix2(qureg, 0, u).unwrap();
    apply_matrix2(qureg, 1, u).unwrap();

    apply_matrix2(qureg, -1, u).unwrap_err();
    apply_matrix2(qureg, 2, u).unwrap_err();
}

#[test]
fn apply_matrix4_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);

    let u = &ComplexMatrix4::new(
        [
            [0., 1., 0., 0.],
            [1., 0., 0., 0.],
            [0., 0., 1., 0.],
            [0., 0., 0., 1.],
        ],
        [
            [0., 0., 0., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 0.],
        ],
    );
    apply_matrix4(qureg, 0, 1, u).unwrap();
    apply_matrix4(qureg, 1, 0, u).unwrap();

    apply_matrix4(qureg, 0, 0, u).unwrap_err();
    apply_matrix4(qureg, 1, 1, u).unwrap_err();

    apply_matrix4(qureg, -1, 1, u).unwrap_err();
    apply_matrix4(qureg, 3, 1, u).unwrap_err();
    apply_matrix4(qureg, 0, 3, u).unwrap_err();
    apply_matrix4(qureg, 0, -3, u).unwrap_err();

    apply_matrix4(qureg, 3, -3, u).unwrap_err();
}

#[test]
fn apply_qft_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();
    init_zero_state(qureg);

    apply_qft(qureg, &[0, 1]).unwrap();
    apply_qft(qureg, &[1, 0]).unwrap();
    apply_qft(qureg, &[1, 2]).unwrap();
    apply_qft(qureg, &[0, 2]).unwrap();
}

#[test]
fn apply_qft_02() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();
    init_zero_state(qureg);

    apply_qft(qureg, &[0, 0]).unwrap_err();
    apply_qft(qureg, &[1, 1]).unwrap_err();
    apply_qft(qureg, &[-1, 0]).unwrap_err();
    apply_qft(qureg, &[4, 0]).unwrap_err();
}

#[test]
fn apply_projector_01() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);
    apply_projector(qureg, 0, 0).unwrap();
    init_zero_state(qureg);
    apply_projector(qureg, 1, 0).unwrap();
    init_zero_state(qureg);
    apply_projector(qureg, 0, 1).unwrap();
    init_zero_state(qureg);
    apply_projector(qureg, 1, 1).unwrap();
}

#[test]
fn apply_projector_02() {
    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(2, env).unwrap();
    init_zero_state(qureg);
    apply_projector(qureg, 0, -1).unwrap_err();
    apply_projector(qureg, 0, 3).unwrap_err();
    apply_projector(qureg, 2, 0).unwrap_err();
    apply_projector(qureg, -1, 0).unwrap_err();
}

#[test]
fn multi_rotate_pauli_01() {
    use PauliOpType::PAULI_X;

    let env = &QuestEnv::new();
    let qureg = &mut Qureg::try_new(3, env).unwrap();
    init_zero_state(qureg);

    multi_rotate_pauli(qureg, &[0, 1], &[PAULI_X, PAULI_X], 0.).unwrap();
    multi_rotate_pauli(qureg, &[1, 2], &[PAULI_X, PAULI_X], 0.).unwrap();
    multi_rotate_pauli(qureg, &[2, 0], &[PAULI_X, PAULI_X], 0.).unwrap();

    multi_rotate_pauli(qureg, &[0, 0], &[PAULI_X, PAULI_X], 0.).unwrap_err();
    multi_rotate_pauli(qureg, &[1, 1], &[PAULI_X, PAULI_X], 0.).unwrap_err();
    multi_rotate_pauli(qureg, &[2, 2], &[PAULI_X, PAULI_X], 0.).unwrap_err();

    multi_rotate_pauli(qureg, &[0, 3], &[PAULI_X, PAULI_X], 0.).unwrap_err();
    multi_rotate_pauli(qureg, &[-1, 0], &[PAULI_X, PAULI_X], 0.).unwrap_err();
    multi_rotate_pauli(qureg, &[0, 1, 2], &[PAULI_X, PAULI_X], 0.).unwrap_err();
    multi_rotate_pauli(qureg, &[0, 1, 2, 3], &[PAULI_X, PAULI_X], 0.)
        .unwrap_err();
    multi_rotate_pauli(qureg, &[0, 1], &[PAULI_X], 0.).unwrap_err();
}
