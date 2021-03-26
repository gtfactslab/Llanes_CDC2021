#include "asif_util.hpp"
#include <chrono>
#include <iostream>

using namespace asif;
using namespace matrix;

ASIF::ASIF(OSQPWorkspace *workspace, float gravity, float mass, float moment_inertia_x, float spring_saturation,
           float spring_constant, float asif_alpha, float k_angle_error, float k_roll_rate_damp, float k_velocity_damp,
           float safe_distance_y, float safe_distance_z, float quad_max_thrust, float quad_lateral_moment_arm,
           float quad_longitudinal_moment_arm, float rotor_moment_constant)
           : osqp_workspace_(workspace), gravity_(gravity), mass_(mass), moment_inertia_x_(moment_inertia_x),
             spring_saturation_(spring_saturation), spring_constant_(spring_constant), asif_alpha_(asif_alpha),k_angle_error_(k_angle_error),
             k_roll_rate_damp_(k_roll_rate_damp), k_velocity_damp_(k_velocity_damp), safe_distance_y_(safe_distance_y),
             safe_distance_z_(safe_distance_z), quad_max_thrust_(quad_max_thrust), quad_lateral_moment_arm_(quad_lateral_moment_arm),
             quad_longitudinal_moment_arm_(quad_longitudinal_moment_arm), rotor_moment_constant_(rotor_moment_constant)
{
    float xbar[10] = {0.F,0.F,0.F,0.F,0.F,0.F,0.F,0.F,safe_distance_y_,safe_distance_z_};
    (void)memcpy(xbar_,xbar,sizeof(xbar));

    float MatrixP[10][10] = {{10.69638864,-9.530123057,1.983401036e-26,-2.162431967e-26,3.365348361,-2.059346592,0.06690696098,-0.05997243717,-4.324192895,-8.8975e-27},
                             {-9.530123057,10.69638864,-1.982571506e-26,2.161537523e-26,-2.059346592,3.365348361,-0.05997243717,0.06690696098,4.324192895,8.8952e-27},
                             {1.983401036e-26,-1.982571506e-26,10.57589286,-9.472321429,4.279724229e-27,-4.301816591e-27,1.246211258e-28,-1.246136893e-28,8.61728179e-27,-4.2917},
                             {-2.162431967e-26,2.161537523e-26,-9.472321429,10.57589286,-4.667164626e-27,4.690079732e-27,-1.357589719e-28,1.357652084e-28,-9.067626648e-27,4.2917},
                             {3.365348361,-2.059346592,4.279724229e-27,-4.667164626e-27,23.78350289,-0.4443626285,0.02146725789,-0.01294033862,-1.070215977,-1.9233e-27},
                             {-2.059346592,3.365348361,-4.301816591e-27,4.690079732e-27,-0.4443626285,23.78350289,-0.01294033862,0.02146725789,1.070215977,1.9396e-27},
                             {0.06690696098,-0.05997243717,1.246211258e-28,-1.357589719e-28,0.02146725789,-0.01294033862,0.01518773389,-0.0003768873623,-0.02696759259,-5.5575e-29},
                             {-0.05997243717,0.06690696098,-1.246136893e-28,1.357652084e-28,-0.01294033862,0.02146725789,-0.0003768873623,0.01518773389,0.02696759259,5.558e-29},
                             {-4.324192895,4.324192895,8.61728179e-27,-9.067626648e-27,-1.070215977,1.070215977,-0.02696759259,0.02696759259,3.120665475,-1.5939e-27},
                             {-8.897478474e-27,8.895150698e-27,-4.291666667,4.291666667,-1.923329526e-27,1.939595904e-27,-5.557529102e-29,5.558038475e-29,-1.593870076e-27,3.1123}};

    (void)memcpy(Pf_, MatrixP, 10*10*sizeof(MatrixP[0][0]));
    (void)memcpy(&(matrix_Pf_(0,0)), MatrixP, 10*10*sizeof(MatrixP[0][0]));

    matrix_Sf_.setZero();
    matrix_Sf_(8,8) = 1;
    matrix_Sf_(9,9) = 1;

    getCorners(w_min_, w_max_, w_corners_);
    phi_ = new embed_state_t[backup_steps_];
    getDerivativeCorners();
}
ASIF::~ASIF()
{
    delete [] phi_;
}

int ASIF::QP(const quad_states_s *quads, control_s *controls, Vector<float,20> &phi_max, float &psi, int &phi_max_time, Vector<float,4> &force_moments)
{
    Vector<float,10> f_x;
    Matrix<float,10,4> g_x;

    f_x.setZero();
    f_x(2) = gravity_;
    f_x(3) = gravity_;
    f_x(4) = quads[0].wx;
    f_x(5) = quads[1].wx;
    f_x(8) = quads[1].vy - quads[0].vy;
    f_x(9) = quads[1].vz - quads[0].vz;

    g_x.setZero();

    g_x(0,0) = sinf(quads[0].phi)/mass_;
    g_x(1,1) = sinf(quads[1].phi)/mass_;
    g_x(2,0) = -cosf(quads[0].phi)/mass_;
    g_x(3,1) = -cosf(quads[1].phi)/mass_;
    g_x(6,2) = 1/moment_inertia_x_;
    g_x(7,3) = 1/moment_inertia_x_;

    states_t x;
    x[0] = quads[0].vy;
    x[1] = quads[1].vy;
    x[2] = quads[0].vz;
    x[3] = quads[1].vz;
    x[4] = quads[0].phi;
    x[5] = quads[1].phi;
    x[6] = quads[0].wx;
    x[7] = quads[1].wx;
    x[8] = quads[1].y - quads[0].y;
    x[9] = quads[1].z - quads[0].z;

    float hs_min_backup[backup_barrier_steps_];
    auto t1 = std::chrono::high_resolution_clock::now();
    phiE(x, x); // update private class member phi_ : embedding trajectory
    auto t2 = std::chrono::high_resolution_clock::now();

    int backup_times_backup[backup_barrier_steps_];
    if (QP_initialized && max_time_ind_backup_ != -1) {
        backupTrajectoryDistribution(max_time_ind_backup_, backup_steps_, backup_times_backup);
        for (int i = 0; i < backup_barrier_steps_; i++) {
        }
    } else {
        for (int i = 0; i < backup_barrier_steps_; i++) {
            backup_times_backup[i] = i*(backup_steps_ - 1) / (backup_barrier_steps_ - 1);
        }
    }

    auto t3 = std::chrono::high_resolution_clock::now();
    // Backup set max of min
    max_time_ind_backup_ = 0;
    float Psi_backup = -INFINITY;
    for (int i = 0; i < backup_barrier_steps_; i++) {
        hSamSoftMin(phi_[backup_times_backup[i]],&(hs_min_backup[i]));
        if(hs_min_backup[i] > Psi_backup){
            Psi_backup = hs_min_backup[i];
            max_time_ind_backup_ = backup_times_backup[i];
        }
    }
    auto t4 = std::chrono::high_resolution_clock::now();
    phi_max_time = max_time_ind_backup_;
    phi_max = phi_[max_time_ind_backup_];
    psi = Psi_backup;

    // Safe set min of min
    float Psi_unsafe = INFINITY;
    if ((QP_initialized && min_time_ind_unsafe_ != -1) || (max_time_ind_backup_ <= backup_barrier_steps_)) {
        min_time_ind_unsafe_ = 0;
        int backup_times_unsafe[max_time_ind_backup_+1];
        float hs_min_unsafe[max_time_ind_backup_+1];
        for (int i = 0; i < max_time_ind_backup_+1; i++) {
            backup_times_unsafe[i] = i;
        }

        for (int i = 0; i < max_time_ind_backup_+1; i++) {
            hSamSoftMinUnsafe(phi_[backup_times_unsafe[i]],&(hs_min_unsafe[i]));
            if(hs_min_unsafe[i] < Psi_unsafe){
                Psi_unsafe = hs_min_unsafe[i];
                min_time_ind_unsafe_ = backup_times_unsafe[i];
            }
        }
    } else {
        min_time_ind_unsafe_ = 0;
        int backup_times_unsafe[backup_barrier_steps_];
        float hs_min_unsafe[backup_barrier_steps_];
        if (min_time_ind_unsafe_ >= max_time_ind_backup_) {
            backupTrajectoryDistribution(max_time_ind_backup_ , max_time_ind_backup_, backup_times_unsafe);
        } else {
            backupTrajectoryDistribution(min_time_ind_unsafe_ , max_time_ind_backup_, backup_times_unsafe);
        }

        for (int i = 0; i < backup_barrier_steps_; i++) {
            hSamSoftMinUnsafe(phi_[backup_times_unsafe[i]],&(hs_min_unsafe[i]));

            if(hs_min_unsafe[i] < Psi_unsafe){
                Psi_unsafe = hs_min_unsafe[i];
                min_time_ind_unsafe_ = backup_times_unsafe[i];
            }
        }
    }
    auto t5 = std::chrono::high_resolution_clock::now();
    Matrix<float,20,10> QMatrix_backup;
    Matrix<float,20,10> QMatrix_unsafe;


    if (max_time_ind_backup_ > min_time_ind_unsafe_) {
        Matrix<float,10,10> Q1;
        Q1.setIdentity();
        Matrix<float,10,10> Q2;
        Q2.setIdentity();

        for (int k = 0; k < min_time_ind_unsafe_; k++) {
            states_t phi_x;
            states_t phi_xh;

            (void)memcpy(phi_x, &(phi_[k+1](0)),10*sizeof(phi_[k](0)));
            (void)memcpy(phi_xh, &(phi_[k+1](10)),10*sizeof(phi_[k](0)));

            Q1 = jacobian(phi_x)*Q1*dt_backup_ + Q1;
            Q2 = jacobian(phi_xh)*Q2*dt_backup_ + Q2;
        }

        (void)memcpy(&(QMatrix_unsafe(0,0)), &(Q1(0,0)),10*10*sizeof(Q1(0,0)));
        (void)memcpy(&(QMatrix_unsafe(10,0)), &(Q2(0,0)),10*10*sizeof(Q2(0,0)));
        auto t10 = std::chrono::high_resolution_clock::now();
        auto duration_unsafe_Q = std::chrono::duration_cast<std::chrono::microseconds>(t10 - t5);
        std::cout << duration_unsafe_Q.count() << " Q_unsafe\n";

        for (int k = min_time_ind_unsafe_; k < max_time_ind_backup_; k++) {
            states_t phi_x;
            states_t phi_xh;

            (void)memcpy(phi_x, &(phi_[k+1](0)),10*sizeof(phi_[k](0)));
            (void)memcpy(phi_xh, &(phi_[k+1](10)),10*sizeof(phi_[k](0)));

            Q1 = jacobian(phi_x)*Q1*dt_backup_ + Q1;
            Q2 = jacobian(phi_xh)*Q2*dt_backup_ + Q2;
        }

        (void)memcpy(&(QMatrix_backup(0,0)), &(Q1(0,0)),10*10*sizeof(Q1(0,0)));
        (void)memcpy(&(QMatrix_backup(10,0)), &(Q2(0,0)),10*10*sizeof(Q2(0,0)));

    } else {
        Matrix<float,10,10> Q1;
        Q1.setIdentity();
        Matrix<float,10,10> Q2;
        Q2.setIdentity();

        for (int k = 0; k < max_time_ind_backup_; k++) {
            states_t phi_x;
            states_t phi_xh;

            (void)memcpy(phi_x, &(phi_[k+1](0)),10*sizeof(phi_[k](0)));
            (void)memcpy(phi_xh, &(phi_[k+1](10)),10*sizeof(phi_[k](0)));

            Q1 = jacobian(phi_x)*Q1*dt_backup_ + Q1;
            Q2 = jacobian(phi_xh)*Q2*dt_backup_ + Q2;
        }
        (void)memcpy(&(QMatrix_unsafe(0,0)), &(Q1(0,0)),10*10*sizeof(Q1(0,0)));
        (void)memcpy(&(QMatrix_unsafe(10,0)), &(Q2(0,0)),10*10*sizeof(Q2(0,0)));
        (void)memcpy(&(QMatrix_backup(0,0)), &(Q1(0,0)),10*10*sizeof(Q1(0,0)));
        (void)memcpy(&(QMatrix_backup(10,0)), &(Q2(0,0)),10*10*sizeof(Q2(0,0)));
    }

    auto t12 = std::chrono::high_resolution_clock::now();
    Matrix<float,1,10> DPsi_backup = hSamSoftMinGrad(phi_[max_time_ind_backup_])*QMatrix_backup;
    auto t11 = std::chrono::high_resolution_clock::now();
    Matrix<float,1,10> DPsi_unsafe = hSamSoftMinGradUnsafe(phi_[min_time_ind_unsafe_])*QMatrix_unsafe;
    auto t6 = std::chrono::high_resolution_clock::now();

    float Tau1_cmd = controls[0].f1 + controls[0].f2 + controls[0].f3 + controls[0].f4;
    float Tau2_cmd = controls[1].f1 + controls[1].f2 + controls[1].f3 + controls[1].f4;
    float Mx1_cmd = quad_lateral_moment_arm_*(-controls[0].f1 + controls[0].f2 + controls[0].f3 - controls[0].f4);
    float Mx2_cmd = quad_lateral_moment_arm_*(-controls[1].f1 + controls[1].f2 + controls[1].f3 - controls[1].f4);

    float q_new[4] = {-Tau1_cmd, -Tau2_cmd, -Mx1_cmd, -Mx2_cmd};
    float warm_x[4] = {Tau1_cmd, Tau2_cmd, Mx1_cmd, Mx2_cmd};
    float ub_new[132];
    float A_new[520];

    Matrix<float,1,4> Ax_backup(-DPsi_backup*g_x);
    Matrix<float,1,4> Ax_unsafe(-DPsi_unsafe*g_x);

    int odd_idx = 0;
    for (int k = 0; k < 64; k++) {
        Vector<float,10> w_temp;
        w_temp(0) = w_corners_[0][k];
        w_temp(1) = w_corners_[1][k];
        w_temp(2) = w_corners_[2][k];
        w_temp(3) = w_corners_[3][k];
        w_temp(6) = w_corners_[4][k];
        w_temp(7) = w_corners_[5][k];

        odd_idx = 2*k;
        ub_new[odd_idx] = (DPsi_backup*(f_x + w_temp) + asif_alpha_*copysign(1.0,Psi_backup)*powf(abs(Psi_backup),2.0))(0,0);
        ub_new[odd_idx+1] = (DPsi_unsafe*(f_x + w_temp) + asif_alpha_*copysign(1.0,Psi_unsafe)*powf(abs(Psi_unsafe),1.0))(0,0);
    }

    ub_new[128] = quad_max_thrust_;
    ub_new[129] = quad_max_thrust_;
    ub_new[130] = quad_max_thrust_;
    ub_new[131] = quad_max_thrust_;


    for (int Axindx = 0; Axindx < 520; Axindx+=2) {
        A_new[Axindx] = Ax_backup(0,Axindx / 130);
        A_new[Axindx+1] = Ax_unsafe(0,Axindx / 130);
    }

    A_new[128] = 0.25F;
    A_new[129] = 0.25F;

    A_new[258] = 0.25F;
    A_new[259] = 0.25F;

    A_new[388] = -0.25F/quad_lateral_moment_arm_;
    A_new[389] =  0.25F/quad_lateral_moment_arm_;

    A_new[518] = -0.25F/quad_lateral_moment_arm_;
    A_new[519] =  0.25F/quad_lateral_moment_arm_;

    if (osqp_update_A(osqp_workspace_, A_new, OSQP_NULL, osqp_workspace_->data->A->nzmax)) {
        return 0;
    }
    if (osqp_update_upper_bound(osqp_workspace_, ub_new)) {
        return 0;
    }
    if (osqp_update_lin_cost(osqp_workspace_, q_new)){
        return 0;
    }

    osqp_warm_start_x(osqp_workspace_, warm_x);

    if (osqp_solve(osqp_workspace_)){
        return 0;
    }

    auto t7 = std::chrono::high_resolution_clock::now();

    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2);
    auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3);
    auto duration4 = std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4);
    auto duration5 = std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5);
    auto duration6 = std::chrono::duration_cast<std::chrono::microseconds>(t7 - t6);
    auto duration11 = std::chrono::duration_cast<std::chrono::microseconds>(t6 - t11);
    auto duration12 = std::chrono::duration_cast<std::chrono::microseconds>(t11 - t12);

    std::cout << duration1.count() << "\n";
    std::cout << duration2.count() << "\n";
    std::cout << duration3.count() << "\n";
    std::cout << duration4.count() << "\n";
    std::cout << duration5.count() << "\n";
    std::cout << duration6.count() << "\n";
    std::cout << duration11.count() << " DPsi_unsafe\n";
    std::cout << duration12.count() << " DPsi_backup\n";
    force_moments(0) = osqp_workspace_->solution->x[0];
    force_moments(1) = osqp_workspace_->solution->x[1];
    force_moments(2) = osqp_workspace_->solution->x[2];
    force_moments(3) = osqp_workspace_->solution->x[3];

    float Tau1 = osqp_workspace_->solution->x[0] / cos(quads[0].theta);
    float Tau2 = osqp_workspace_->solution->x[1] / cos(quads[1].theta);
    float Mx1 = osqp_workspace_->solution->x[2];
    float Mx2 = osqp_workspace_->solution->x[3];
    float My1 = -(quads[0].theta) - k_roll_rate_damp_ * quads[0].wy;
    float My2 = -(quads[1].theta) - k_roll_rate_damp_ * quads[1].wy;
    float Mz1 = -(quads[0].psi) - k_roll_rate_damp_ * quads[0].wz;
    float Mz2 = -(quads[1].psi) - k_roll_rate_damp_ * quads[1].wz;

    controls[0].f1 = 0.25F*(-Mx1/quad_lateral_moment_arm_ + My1/quad_longitudinal_moment_arm_ + Mz1/rotor_moment_constant_ + Tau1);
    controls[0].f2 = 0.25F*(Mx1/quad_lateral_moment_arm_ - My1/quad_longitudinal_moment_arm_ + Mz1/rotor_moment_constant_ + Tau1);
    controls[0].f3 = 0.25F*(Mx1/quad_lateral_moment_arm_ + My1/quad_longitudinal_moment_arm_ - Mz1/rotor_moment_constant_ + Tau1);
    controls[0].f4 = 0.25F*(-Mx1/quad_lateral_moment_arm_ - My1/quad_longitudinal_moment_arm_ - Mz1/rotor_moment_constant_ + Tau1);

    controls[1].f1 = 0.25F*(-Mx2/quad_lateral_moment_arm_ + My2/quad_longitudinal_moment_arm_ + Mz2/rotor_moment_constant_ + Tau2);
    controls[1].f2 = 0.25F*(Mx2/quad_lateral_moment_arm_ - My2/quad_longitudinal_moment_arm_ + Mz2/rotor_moment_constant_ + Tau2);
    controls[1].f3 = 0.25F*(Mx2/quad_lateral_moment_arm_ + My2/quad_longitudinal_moment_arm_ - Mz2/rotor_moment_constant_ + Tau2);
    controls[1].f4 = 0.25F*(-Mx2/quad_lateral_moment_arm_ - My2/quad_longitudinal_moment_arm_ - Mz2/rotor_moment_constant_ + Tau2);

    printf("Status:                %s\n", (osqp_workspace_)->info->status);
    printf("Number of iterations:  %d\n", (int)((osqp_workspace_)->info->iter));
    printf("obj_val:               %f\n", (osqp_workspace_)->info->obj_val);
    printf("prim_res:              %f\n", (osqp_workspace_)->info->pri_res);

    printf("Tau1_cmd: %f, Tau2_cmd: %f, Mx1_cmd: %f, Mx2_cmd: %f \n", Tau1_cmd,
                                                                      Tau2_cmd,
                                                                      Mx1_cmd,
                                                                      Mx2_cmd);
    printf("Tau1: %f, Tau2: %f, Mx1: %f, Mx2: %f \n", osqp_workspace_->solution->x[0],
                                                      osqp_workspace_->solution->x[1],
                                                      osqp_workspace_->solution->x[2],
                                                      osqp_workspace_->solution->x[3]);
    (void)printf("Psi_backup:                          %.10f\n", Psi_backup);
    (void)printf("Psi_safe:                          %.10f\n", Psi_unsafe);


    int solution_solved = (OSQP_SOLVED == osqp_workspace_->info->status_val);

    if (solution_solved) {
        QP_initialized = true;
    }

    return solution_solved;
}

void ASIF::hSamSoftMin(embed_state_t phi, float *hs_min)
{
    float m = -100000000.F;
    state_space_corners_t state_corners;
    states_t x;
    states_t xhat;
    (void)memcpy(x,&(phi(0)),10*sizeof(phi(0)));
    (void)memcpy(xhat,&(phi(10)),10*sizeof(phi(0)));

    getCorners(x, xhat, state_corners);

    float hs[1024];
    hsamCorners(state_corners, hs);
    int max_indx = 0;
    for (int i = 0; i < 1024; i++) {
        if(hs[i]*m > hs[max_indx]*m){
            max_indx = i;
        }
    }
    float xstar = hs[max_indx]*m;
    float sum = 0.F;

    for (float h : hs) {
        sum += expf(h*m - xstar);
    }
    *hs_min = 1/m*(xstar + log(sum));
}

void ASIF::hSamSoftMinUnsafe(embed_state_t phi, float *hs_min) {
    float m = -100000000.F;
    alpha_beta_corners_t state_corners;
    states_t x;
    states_t xhat;
    (void)memcpy(x,&(phi(0)),10*sizeof(phi(0)));
    (void)memcpy(xhat,&(phi(10)),10*sizeof(phi(0)));

    getCornersUnsafe(x, xhat, state_corners);

    float hs[4];
    hsamCornersUnsafe(state_corners, hs);
    int max_indx = 0;
    for (int i = 0; i < 4; i++) {
        if(hs[i]*m > hs[max_indx]*m){
            max_indx = i;
        }
    }
    float xstar = hs[max_indx]*m;
    float sum = 0.F;

    for (float h : hs) {
        sum += expf(h*m - xstar);
    }
    *hs_min = 1/m*(xstar + log(sum));
}

Matrix<float,1,20> ASIF::hSamSoftMinGrad(embed_state_t phi)
{
    float m = -100000000.F;
    state_space_corners_t state_corners;
    states_t x;
    states_t xhat;
    Matrix<float,1,20> out;
    Matrix<float,10,1> xbar(xbar_);
    (void)memcpy(x,&(phi(0)),10*sizeof(phi(0)));
    (void)memcpy(xhat,&(phi(10)),10*sizeof(phi(0)));
    getCorners(x, xhat, state_corners);

    Matrix<float,10,1024> corners(state_corners);

    float hs[1024];
    hsamCorners(state_corners, hs);
    int max_indx = 0;
    for (int i = 0; i < 1024; i++) {
        if(hs[i]*m > hs[max_indx]*m){
            max_indx = i;
        }
    }
    float xstar = hs[max_indx]*m;
    float sum = 0.F;

    for (float h : hs) {
        sum += expf(h*m - xstar);
    }
    for (int j = 0; j < 1024; j++) {
        Matrix<float,10,1> corner(corners.slice<10,1>(0,j));
        out = out - 2 * expf(m*hs[j]-xstar) * (corner - xbar).transpose() * matrix_Pf_ * derivative_corners_[j];
    }
    out = 1/sum*out;
    return out;
}

Matrix<float,1,20> ASIF::hSamSoftMinGradUnsafe(embed_state_t phi) {
    float m = -100000000.F;
    alpha_beta_corners_t alpha_beta_corners;
    states_t x;
    states_t xhat;
    Matrix<float,1,20> out;
    Matrix<float,10,1> xbar(xbar_);
    (void)memcpy(x,&(phi(0)),10*sizeof(phi(0)));
    (void)memcpy(xhat,&(phi(10)),10*sizeof(phi(0)));
    getCornersUnsafe(x, xhat, alpha_beta_corners);

    Matrix<float,10,4> corners(alpha_beta_corners);

    float hs[4];
    hsamCornersUnsafe(alpha_beta_corners, hs);
    int max_indx = 0;
    for (int i = 0; i < 4; i++) {
        if(hs[i]*m > hs[max_indx]*m){
            max_indx = i;
        }
    }
    float xstar = hs[max_indx]*m;
    float sum = 0.F;

    for (float h : hs) {
        sum += expf(h*m - xstar);
    }
    for (int j = 0; j < 4; j++) {
        Matrix<float,10,1> corner(corners.slice<10,1>(0,j));
        out = out - 2 * expf(m*hs[j]-xstar) * (corner - xbar).transpose() * matrix_Sf_ * derivative_corners_unsafe_[j];
    }
    out = 1/sum*out;
    return out;
}

void ASIF::getCorners(const float *w_min, const float *w_max, disturbance_corners_t corners)
{
    for(int j = 0; j < 64; j++) {
        for(int i = 0; i < 6; i++) {
            if ((j & (0b000001 << i)) == (0b000001 << i)) {
                corners[i][j] = w_max[i];
            } else {
                corners[i][j] = w_min[i];
            }
        }
    }
}

void ASIF::getCorners(states_t x_min, states_t x_max, state_space_corners_t corners)
{
    float corners_[10][1024] = {
    				  { x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],\
                        x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0],x_min[0],x_max[0]},\
                      { x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1],\
                        x_min[1],x_min[1],x_max[1],x_max[1],x_min[1],x_min[1],x_max[1],x_max[1]},\
                      { x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2],\
                        x_min[2],x_min[2],x_min[2],x_min[2],x_max[2],x_max[2],x_max[2],x_max[2]},\
                      { x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],\
                        x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],x_min[3],\
                        x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3],x_max[3]},\
                      { x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],x_min[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],\
                        x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4],x_max[4]},\
                      { x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],x_min[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],\
                        x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5],x_max[5]},\
                      { x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],x_min[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],\
                        x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6],x_max[6]},\
                      { x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],x_min[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],\
                        x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7],x_max[7]},\
                      { x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],x_min[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],\
                        x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8],x_max[8]},\
                      { x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],x_min[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],\
                        x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9],x_max[9]}};

        memcpy(&corners[0][0],&corners_[0][0],1024*10*sizeof(corners_[0][0]));
}

void ASIF::getCornersUnsafe(states_t x_min, states_t x_max, alpha_beta_corners_t corners) {
    float corners_[10][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},\
                             {0,0,0,0},{x_min[8],x_max[8],x_min[8],x_max[8]},{x_min[9],x_min[9],x_max[9],x_max[9]}};
    memcpy(&corners[0][0],&corners_[0][0],4*10*sizeof(corners_[0][0]));
}

void ASIF::hsamCorners(state_space_corners_t corners, float *hs)
{
    for (int k = 0; k < 1024; k++){
        hs[k] = barrier_constant_backup_;
        for (int i = 0; i < 10; i++) {
            hs[k] -= (corners[i][k] - xbar_[i])*((corners[0][k] - xbar_[0])*Pf_[i][0] + (corners[1][k] - xbar_[1])*Pf_[i][1]\
                    + (corners[2][k] - xbar_[2])*Pf_[i][2] + (corners[3][k] - xbar_[3])*Pf_[i][3]\
                    + (corners[4][k] - xbar_[4])*Pf_[i][4]  + (corners[5][k] - xbar_[5])*Pf_[i][5]\
                    + (corners[6][k] - xbar_[6])*Pf_[i][6]  + (corners[7][k] - xbar_[7])*Pf_[i][7]\
                    + (corners[8][k] - xbar_[8])*Pf_[i][8]  + (corners[9][k] - xbar_[9])*Pf_[i][9]);
        }
    }
}

void ASIF::hsamCornersUnsafe(alpha_beta_corners_t corners, float *hs)
{
    for (int k = 0; k < 4; k++){
        hs[k] = powf(safe_radius_,2.0F) - powf(corners[8][k] - xbar_[8],2.0F) - powf(corners[9][k] - xbar_[9],2.0F);
    }
}

Matrix<float,10,10> ASIF::jacobian(states_t x)
{
    Matrix<float,10,10> J;
    float alpha_tanh = tanhf(spring_saturation_*(x[8] - safe_distance_y_));
    float beta_tanh = tanhf(spring_saturation_*(x[9] - safe_distance_z_));
    float fy1 = spring_constant_ * alpha_tanh - k_velocity_damp_ * x[0];
    float fy2 = -spring_constant_ * alpha_tanh - k_velocity_damp_ * x[1];
    float fz1 = spring_constant_ * beta_tanh - mass_ * gravity_ - k_velocity_damp_ * x[2];
    float fz2 = -spring_constant_ * beta_tanh - mass_ * gravity_ - k_velocity_damp_ * x[3];

    J(0,0) = -k_velocity_damp_/mass_*powf(sinf(x[4]),2.F);
    J(0,1) = 0.F;
    J(0,2) = k_velocity_damp_/mass_*sinf(x[4])*cosf(x[4]);
    J(0,3) = 0.F;
    J(0,4) = 1.F/mass_*(fy1 * sinf(2 * x[4]) - fz1 * cosf(2.F * x[4]));
    J(0,5) = 0.F;
    J(0,6) = 0.F;
    J(0,7) = 0.F;
    J(0,8) = spring_constant_*spring_saturation_/mass_*powf(sinf(x[4]),2.F)*(1.F-powf(alpha_tanh,2.F));
    J(0,9) = -spring_constant_*spring_saturation_/mass_*sinf(x[4])*cosf(x[4])*(1.F-powf(beta_tanh,2.F));

    J(1,0) = 0.F;
    J(1,1) = -k_velocity_damp_/mass_*powf(sinf(x[5]),2.F);
    J(1,2) = 0.F;
    J(1,3) = k_velocity_damp_/mass_*sinf(x[5])*cosf(x[5]);
    J(1,4) = 0.F;
    J(1,5) = 1.F/mass_*(fy2 * sinf(2.F * x[5]) - fz2 * cosf(2.F * x[5]));
    J(1,6) = 0.F;
    J(1,7) = 0.F;
    J(1,8) = -spring_constant_*spring_saturation_/mass_*powf(sinf(x[5]),2)*(1.F-powf(alpha_tanh,2.F));
    J(1,9) = spring_constant_*spring_saturation_/mass_*sinf(x[5])*cosf(x[5])*(1.F-powf(beta_tanh,2.F));

    J(2,0) = k_velocity_damp_/mass_*sinf(x[4])*cosf(x[4]);
    J(2,1) = 0.F;
    J(2,2) = -k_velocity_damp_/mass_*powf(cosf(x[4]),2.F);
    J(2,3) = 0.F;
    J(2,4) = -1.F/mass_*(fy1 * cosf(2.F * x[4]) + fz1 * sinf(2.F * x[4]));
    J(2,5) = 0.F;
    J(2,6) = 0.F;
    J(2,7) = 0.F;
    J(2,8) = -spring_constant_*spring_saturation_/mass_*sinf(x[4])*cosf(x[4])*(1-powf(alpha_tanh,2.F));
    J(2,9) = spring_constant_*spring_saturation_/mass_*powf(cosf(x[4]),2)*(1-powf(beta_tanh,2.F));

    J(3,0) = 0.F;
    J(3,1) = k_velocity_damp_/mass_*sinf(x[5])*cosf(x[5]);
    J(3,2) = 0.F;
    J(3,3) = -k_velocity_damp_/mass_*powf(cosf(x[5]),2.F);
    J(3,4) = 0.F;
    J(3,5) = -1.F/mass_*(fy2 * cosf(2 * x[5]) + fz2 * sinf(2 * x[5]));
    J(3,6) = 0.F;
    J(3,7) = 0.F;
    J(3,8) = spring_constant_*spring_saturation_/mass_*sinf(x[5])*cosf(x[5])*(1-powf(alpha_tanh,2.F));
    J(3,9) = -spring_constant_*spring_saturation_/mass_*powf(cosf(x[5]),2)*(1-powf(beta_tanh,2.F));

    J(4,0) = 0.F;
    J(4,1) = 0.F;
    J(4,2) = 0.F;
    J(4,3) = 0.F;
    J(4,4) = 0.F;
    J(4,5) = 0.F;
    J(4,6) = 1.F;
    J(4,7) = 0.F;
    J(4,8) = 0.F;
    J(4,9) = 0.F;

    J(5,0) = 0.F;
    J(5,1) = 0.F;
    J(5,2) = 0.F;
    J(5,3) = 0.F;
    J(5,4) = 0.F;
    J(5,5) = 0.F;
    J(5,6) = 0.F;
    J(5,7) = 1.F;
    J(5,8) = 0.F;
    J(5,9) = 0.F;

    J(6,0) = -k_velocity_damp_*k_angle_error_/moment_inertia_x_*cosf(x[4]);
    J(6,1) = 0.F;
    J(6,2) = -k_velocity_damp_*k_angle_error_/moment_inertia_x_*sinf(x[4]);
    J(6,3) = 0.F;
    J(6,4) = k_angle_error_/moment_inertia_x_*(fz1 * cosf(x[4]) - fy1 * sinf(x[4]));
    J(6,5) = 0.F;
    J(6,6) = -k_roll_rate_damp_/moment_inertia_x_;
    J(6,7) = 0.F;
    J(6,8) = k_angle_error_*spring_constant_*spring_saturation_/moment_inertia_x_*(1.F-powf(alpha_tanh,2.F))*cosf(x[4]);
    J(6,9) = k_angle_error_*spring_constant_*spring_saturation_/moment_inertia_x_*(1.F-powf(beta_tanh,2.F))*sinf(x[4]);

    J(7,0) = 0.F;
    J(7,1) = -k_velocity_damp_*k_angle_error_/moment_inertia_x_*cosf(x[5]);
    J(7,2) = 0.F;
    J(7,3) = -k_velocity_damp_*k_angle_error_/moment_inertia_x_*sinf(x[5]);
    J(7,4) = 0.F;
    J(7,5) = k_angle_error_/moment_inertia_x_*(fz2 * cosf(x[5]) - fy2 * sinf(x[5]));
    J(7,6) = 0.F;
    J(7,7) = -k_roll_rate_damp_/moment_inertia_x_;
    J(7,8) = -k_angle_error_*spring_constant_*spring_saturation_/moment_inertia_x_*(1.F-powf(alpha_tanh,2.F))*cosf(x[5]);
    J(7,9) = -k_angle_error_*spring_constant_*spring_saturation_/moment_inertia_x_*(1.F-powf(beta_tanh,2.F))*sinf(x[5]);

    J(8,0) = -1.F;
    J(8,1) = 1.F;
    J(8,2) = 0.F;
    J(8,3) = 0.F;
    J(8,4) = 0.F;
    J(8,5) = 0.F;
    J(8,6) = 0.F;
    J(8,7) = 0.F;
    J(8,8) = 0.F;
    J(8,9) = 0.F;

    J(9,0) = 0.F;
    J(9,1) = 0.F;
    J(9,2) = -1.F;
    J(9,3) = 1.F;
    J(9,4) = 0.F;
    J(9,5) = 0.F;
    J(9,6) = 0.F;
    J(9,7) = 0.F;
    J(9,8) = 0.F;
    J(9,9) = 0.F;

    return J;
}

void ASIF::phiE(states_t x0, states_t xh0){

    memcpy(&(phi_[0](0)),x0, 10*sizeof(x0[0]));
    memcpy(&(phi_[0](10)),xh0, 10*sizeof(xh0[0]));

    for (int i = 0; i < (backup_steps_-1);i++) {
        phi_[i+1] = phi_[i] + dt_backup_ * embeddingDynamics(phi_[i]);
    }
}

embed_state_t ASIF::embeddingDynamics(embed_state_t phi)
{

    embed_state_t out;
    states_t x;
    states_t xh;

    phi.slice<10,1>(0,0).copyTo(x);
    phi.slice<10,1>(10,0).copyTo(xh);

    Vector<float,10> out1(decomp(x,xh,w_min_));
    Vector<float,10> out2(decomp(xh,x,w_max_));

    out(0) = out1(0);
    out(1) = out1(1);
    out(2) = out1(2);
    out(3) = out1(3);
    out(4) = out1(4);
    out(5) = out1(5);
    out(6) = out1(6);
    out(7) = out1(7);
    out(8) = out1(8);
    out(9) = out1(9);
    out(10) = out2(0);
    out(11) = out2(1);
    out(12) = out2(2);
    out(13) = out2(3);
    out(14) = out2(4);
    out(15) = out2(5);
    out(16) = out2(6);
    out(17) = out2(7);
    out(18) = out2(8);
    out(19) = out2(9);
    return out;
}

Vector<float,10> ASIF::decomp(const states_t x, const states_t xh, const float w[6]) const
{
    Vector<float,10> out;
    float y1;
    float y2;
    float y3;
    float y4;
    float y9;
    float y10;
    float za;
    float zb;
    float theta_s;
    float theta_ss;
    float J;
    float q;
    float q_under;
    float q_over;

    // ------------------------- D1 ----------------------------------------

    if (sinf(x[4])*cosf(x[4]) >= 0.F) {
        y3 = x[2];
        y10 = xh[9];
    } else {
        y3 = xh[2];
        y10 = x[9];
    }

    za = spring_constant_/mass_*tanhf(spring_saturation_*(x[8] - safe_distance_y_)) - k_velocity_damp_/mass_*x[0];
    zb = -spring_constant_/mass_*tanhf(spring_saturation_*(y10 - safe_distance_z_)) + gravity_ + k_velocity_damp_/mass_*y3;

    theta_s = 0.5F*atan2f(-zb,za);
    theta_ss = 0.5F*atan2f(za,zb);

    J = za*sinf(2.F*theta_ss) + zb*cosf(2.F*theta_ss);

    q = za*powf(sinf(x[4]),2.F) + zb*cosf(x[4])*sinf(x[4]) + J*(x[4] - xh[4]);
    q_under = za*powf(sinf(theta_s),2.F) + zb/2.F*sinf(2.F*theta_s);
    q_over = za*powf(cosf(theta_s),2.F) - zb/2.F*sinf(2.F*theta_s);

    out(0) = fminf(q_over,fmaxf(q_under,q)) + w[0];

    // ------------------------- D2 ----------------------------------------

    if (sinf(x[5])*cosf(x[5]) >= 0.F) {
        y4 = x[3];
        y10 = x[9];
    } else {
        y4 = xh[3];
        y10 = xh[9];
    }

    za = -spring_constant_/mass_*tanhf(spring_saturation_*(xh[8]-safe_distance_y_)) - k_velocity_damp_/mass_*x[1];
    zb = spring_constant_/mass_*tanhf(spring_saturation_*(y10-safe_distance_z_)) + gravity_ + k_velocity_damp_/mass_*y4;

    theta_s = 0.5F*atan2f(-zb,za);
    theta_ss = 0.5F*atan2f(za,zb);
    J = za*sinf(2.F*theta_ss) + zb*cosf(2.F*theta_ss);

    q = za*powf(sinf(x[5]),2.F) + zb*cosf(x[5])*sinf(x[5]) + J*(x[5] - xh[5]);
    q_under = za*powf(sinf(theta_s),2.F) + zb/2.F*sinf(2.F*theta_s);
    q_over = za*powf(cosf(theta_s),2.F) - zb/2.F*sinf(2.F*theta_s);

    out(1) = fminf(q_over, fmaxf(q_under,q)) + w[1];

    // ------------------------- D3 ----------------------------------------

    if (sinf(x[4])*cosf(x[4]) >= 0.F) {
        y1 = x[0];
        y9 = xh[8];
    } else {
        y1 = xh[0];
        y9 = x[8];
    }

    za = spring_constant_/mass_*tanhf(spring_saturation_*(x[9] - safe_distance_z_)) - gravity_ - k_velocity_damp_/mass_*x[2];
    zb = -spring_constant_/mass_*tanhf(spring_saturation_*(y9 - safe_distance_y_)) + k_velocity_damp_/mass_*y1;

    theta_s = 0.5F*atan2f(zb,za);
    theta_ss = 0.5F*atan2f(-za,zb);

    J = zb*cosf(2.F*theta_ss) - za*sinf(2.F*theta_ss);

    q = za*powf(cosf(x[4]),2.F) + zb*cosf(x[4])*sinf(x[4]) + gravity_ + J*(x[4] - xh[4]);
    q_under = za*powf(sinf(theta_s),2.F) - zb/2.F*sinf(2.F*theta_s) + gravity_;
    q_over = za*powf(cosf(theta_s),2.F) + zb/2.F*sinf(2.F*theta_s) + gravity_;

    out(2) = fminf(q_over, fmaxf(q_under,q)) + w[2];

    // ------------------------- D4 ----------------------------------------

    if (sinf(x[5])*cosf(x[5]) >= 0.F) {
        y2 = x[1];
        y9 = x[8];
    } else {
        y2 = xh[1];
        y9 = xh[8];
    }

    za = -spring_constant_/mass_*tanhf(spring_saturation_*(xh[9]-safe_distance_z_)) - gravity_ - k_velocity_damp_/mass_*x[3];
    zb = spring_constant_/mass_*tanhf(spring_saturation_*(y9-safe_distance_y_)) + k_velocity_damp_/mass_*y2;

    theta_s = 0.5F*atan2f(zb,za);
    theta_ss = 0.5F*atan2f(-za,zb);
    J = zb*cosf(2.F*theta_ss) - za*sinf(2.F*theta_ss);

    q = za*powf(cosf(x[5]),2.F) + zb*cosf(x[5])*sinf(x[5]) + gravity_ + J*(x[5] - xh[5]);
    q_under = za*powf(sinf(theta_s),2.F) - zb/2.F*sinf(2.F*theta_s) + gravity_;
    q_over = za*powf(cosf(theta_s),2.F) + zb/2.F*sinf(2.F*theta_s) + gravity_;

    out(3) = fminf(q_over, fmaxf(q_under,q)) + w[3];

    // ------------------------- D5 and D6 ----------------------------------

    out(4) = x[6];
    out(5) = x[7];

    // ------------------------- D7 ----------------------------------------

    if (cosf(x[4]) <= 0.F) {
        y1 = x[0];
        y9 = xh[8];
    } else {
        y1 = xh[0];
        y9 = x[8];
    }

    if (sinf(x[4]) <= 0.F) {
        y3 = x[2];
        y10 = xh[9];
    } else {
        y3 = xh[2];
        y10 = x[9];
    }

    za = k_angle_error_/moment_inertia_x_*(spring_constant_*tanhf(spring_saturation_*(y10-safe_distance_z_))
            - mass_*gravity_ - k_velocity_damp_*y3);
    zb = k_angle_error_/moment_inertia_x_*(spring_constant_*tanhf(spring_saturation_*(y9-safe_distance_y_))
            - k_velocity_damp_*y1);

    theta_s = atan2f(za,zb);
    theta_ss = -atan2f(za,zb);

    J = za*cosf(theta_ss) - zb*sinf(theta_ss);

    q = za*sinf(x[4]) + zb*cosf(x[4]) + J*(x[4] - xh[4]) - k_roll_rate_damp_/moment_inertia_x_*x[6];
    q_over = za*sinf(theta_s) + zb*cosf(theta_s) - k_roll_rate_damp_/moment_inertia_x_*x[6];
    q_under = -za*sin(theta_s) - zb*cosf(theta_s) - k_roll_rate_damp_/moment_inertia_x_*x[6];

    out(6) = fminf(q_over, fmaxf(q_under,q)) + w[4];


    // ------------------------- D8 ----------------------------------------

    if (cosf(x[5]) <= 0.F) {
        y2 = x[1];
        y9 = x[8];
    } else {
        y2 = xh[1];
        y9 = xh[8];
    }

    if (sinf(x[5]) <= 0.F) {
        y4 = x[3];
        y10 = x[9];
    } else {
        y4 = xh[3];
        y10 = xh[9];
    }

    za = -k_angle_error_/moment_inertia_x_*(spring_constant_*tanhf(spring_saturation_*(y10-safe_distance_z_))
                                           + mass_*gravity_ + k_velocity_damp_*y4);
    zb = -k_angle_error_/moment_inertia_x_*(spring_constant_*tanhf(spring_saturation_*(y9-safe_distance_y_))
                                           + k_velocity_damp_*y2);

    theta_s = atan2f(za,zb);
    theta_ss = -atan2f(za,zb);

    J = za*cosf(theta_ss) - zb*sinf(theta_ss);

    q = za*sinf(x[5]) + zb*cosf(x[5]) + J*(x[5] - xh[5]) - k_roll_rate_damp_/moment_inertia_x_*x[7];
    q_over = za*sinf(theta_s) + zb*cosf(theta_s) - k_roll_rate_damp_/moment_inertia_x_*x[7];
    q_under = -za*sin(theta_s) - zb*cosf(theta_s) - k_roll_rate_damp_/moment_inertia_x_*x[7];

    out(7) = fminf(q_over, fmaxf(q_under,q)) + w[5];

    // ------------------------- D5 and D6 ----------------------------------
    out(8) = x[1] - xh[0];
    out(9) = x[3] - xh[2];

    return out;
}

void ASIF::getDerivativeCorners()
{
    for (int i = 0; i < 1024; i++) {
        derivative_corners_[i].setZero();
        int diag_bin_neg = 1023 - i;
        for (int j = 0; j < 10; j++) {
            if ((i & (0b0000000001 << j)) == (0b0000000001 << j)) {
                derivative_corners_[i](j, j) = 0.F;
                derivative_corners_[i](j, j + 10) = 1.F;
            } else {
                derivative_corners_[i](j, j) = 1.F;
                derivative_corners_[i](j, j + 10) = 0.F;
            }
        }
    }

    derivative_corners_unsafe_[0].setZero();
    derivative_corners_unsafe_[0](8,8) = 1;
    derivative_corners_unsafe_[0](9,9) = 1;
    derivative_corners_unsafe_[1].setZero();
    derivative_corners_unsafe_[1](9,9) = 1;
    derivative_corners_unsafe_[1](8,18) = 1;
    derivative_corners_unsafe_[2].setZero();
    derivative_corners_unsafe_[2](8,8) = 1;
    derivative_corners_unsafe_[2](9,19) = 1;
    derivative_corners_unsafe_[3].setZero();
    derivative_corners_unsafe_[3](8,18) = 1;
    derivative_corners_unsafe_[3](9,19) = 1;

}

void ASIF::backupTrajectoryDistribution(int traj_t, int backup_steps, int *out)
{
    float spacing_factor = 0.01F;
    float corner_power = 7.0F; //must be odd

    int n_right = static_cast<int>(round(0.5F*(-pow((((float)traj_t-(float)backup_steps/2.0F)),corner_power)
            / pow((float)backup_steps/2.0F,corner_power)+1)*(float)backup_barrier_steps_)+1);
    int n_left = backup_barrier_steps_ - n_right + 1;

    auto log_right_step = static_cast<float>((exp(((float)backup_steps-(float)traj_t)*spacing_factor) - 1)/((float)n_right-1));
    auto log_left_step = static_cast<float>((exp(((float)traj_t)*spacing_factor) - 1)/((float)n_left-1));

    for (int i = 0; i < n_left - 1; i++) {
        out[i] = static_cast<int>(round(log(1 + log_left_step*(float)i) / spacing_factor));
    }
    for (int j = 0; j < n_right; j++) {
        out[backup_barrier_steps_ - j - 1] = backup_steps - 1 - static_cast<int>(round(log(1 + log_right_step*(float)j) / spacing_factor));
    }
}

void ASIF::backupControl(const quad_states_s states[2], control_s controls[2])
{
    float f_cmd_x1 =  spring_constant_*tanhf(spring_saturation_*(-states[0].x)) - k_velocity_damp_*states[0].vx;
    float f_cmd_x2 =  spring_constant_*tanhf(spring_saturation_*(-states[1].x)) - k_velocity_damp_*states[1].vx;
    float f_cmd_y1 =  spring_constant_*tanhf(spring_saturation_*((states[1].y - states[0].y) - safe_distance_y_)) - k_velocity_damp_*states[0].vy;
    float f_cmd_y2 =  -spring_constant_*tanhf(spring_saturation_*((states[1].y - states[0].y) - safe_distance_y_)) - k_velocity_damp_*states[1].vy;
    float f_cmd_z1 =  spring_constant_*tanhf(spring_saturation_*((states[1].z - states[0].z) - safe_distance_z_)) - mass_*gravity_ - k_velocity_damp_*states[0].vz;
    float f_cmd_z2 =  -spring_constant_*tanhf(spring_saturation_*((states[1].z - states[0].z) - safe_distance_z_)) - mass_*gravity_ - k_velocity_damp_*states[1].vz;

    Vector3f thrust_dir1 {f_cmd_x1, f_cmd_y1, f_cmd_z1};
    Vector3f thrust_dir2 {f_cmd_x2, f_cmd_y2, f_cmd_z2};

    thrust_dir1.normalize();
    thrust_dir2.normalize();

    float psi1_cmd = 0.F;
    float psi2_cmd = 0.F;

    Dcmf iterm_psi_rotation_matrix1(Eulerf(0.F,0.F,psi1_cmd));
    Dcmf iterm_psi_rotation_matrix2(Eulerf(0.F,0.F,psi1_cmd));

    Vector3f iterm_thrust_frame1_1 = iterm_psi_rotation_matrix1*thrust_dir1;
    Vector3f iterm_thrust_frame1_2 = iterm_psi_rotation_matrix2*thrust_dir2;

    float theta1_cmd = atan2f(-iterm_thrust_frame1_1(0),-iterm_thrust_frame1_1(2));
    float theta2_cmd = atan2f(-iterm_thrust_frame1_2(0),-iterm_thrust_frame1_2(2));

    Dcmf iterm_theta_rotation_matrix1(Eulerf(0.F,theta1_cmd,0.F));
    Dcmf iterm_theta_rotation_matrix2(Eulerf(0.F,theta2_cmd,0.F));

    Vector3f iterm_thrust_frame2_1 = iterm_theta_rotation_matrix1*iterm_thrust_frame1_1;
    Vector3f iterm_thrust_frame2_2 = iterm_theta_rotation_matrix2*iterm_thrust_frame1_2;

    float phi1_cmd = atan2f(iterm_thrust_frame2_1(1),-iterm_thrust_frame2_1(2));
    float phi2_cmd = atan2f(iterm_thrust_frame2_2(1),-iterm_thrust_frame2_2(2));

    float Tau1 = - f_cmd_x1 * cosf(states[0].phi) * sinf(states[0].theta) + f_cmd_y1 * sinf(states[0].phi)
                 - f_cmd_z1 * cosf(states[0].phi) * cosf(states[0].theta);
    float Tau2 = - f_cmd_x2 * cosf(states[1].phi) * sinf(states[1].theta) + f_cmd_y2 * sinf(states[1].phi)
                 - f_cmd_z2 * cosf(states[1].phi) * cosf(states[1].theta);

//    float Mx1 = -k_angle_error_ * (states[0].phi - phi1_cmd) - k_roll_rate_damp_ * states[0].wx;
//    float Mx2 = -k_angle_error_ * (states[1].phi - phi2_cmd) - k_roll_rate_damp_ * states[1].wx;

    float Mx1 = k_angle_error_*f_cmd_z1*sinf(states[0].phi) + k_angle_error_*f_cmd_y1*cosf(states[0].phi) - k_roll_rate_damp_*states[0].wx;
    float Mx2 = k_angle_error_*f_cmd_z2*sinf(states[1].phi) + k_angle_error_*f_cmd_y2*cosf(states[1].phi) - k_roll_rate_damp_*states[1].wx;

    float My1 = -(states[0].theta - theta1_cmd) - k_roll_rate_damp_ * states[0].wy;
    float My2 = -(states[1].theta - theta2_cmd) - k_roll_rate_damp_ * states[1].wy;

    float Mz1 = -(states[0].psi - psi1_cmd) - k_roll_rate_damp_ * states[0].wz;
    float Mz2 = -(states[1].psi - psi2_cmd) - k_roll_rate_damp_ * states[1].wz;

    printf("fcmd_x2: %f,fcmd_y2: %f,fcmd_z2: %f\n", f_cmd_x2, f_cmd_y2, f_cmd_z2);
    printf("phi1: %f, theta1: %f, psi1: %f\n", states[0].phi, states[0].theta, states[0].psi);
    printf("phi1_cmd: %f, theta1_cmd: %f, psi1_cmd: %f\n", phi1_cmd, theta1_cmd, psi1_cmd);
    printf("phi2: %f, theta2: %f, psi2: %f\n", states[1].phi, states[1].theta, states[1].psi);
    printf("phi2_cmd: %f, theta2_cmd: %f, psi2_cmd: %f\n", phi2_cmd, theta2_cmd, psi2_cmd);
    printf("vx2: %f, vy2: %f, zz: %f\n",states[1].vx,states[1].vy,states[1].vz);

    printf("Mx1: %f, My1: %f, Mz1: %f\n", Mx1, My1, Mz1);
    printf("Mx2: %f, My2: %f, Mz2: %f\n", Mx2, My2, Mz2);

    printf("x1: %f, y1: %f, z1: %f\n",states[0].x,states[0].y,states[0].z);
    printf("x2: %f, y2: %f, zz: %f\n\n",states[1].x,states[1].y,states[1].z);

    controls[0].f1 = 0.25F*(-Mx1/quad_lateral_moment_arm_ + My1/quad_longitudinal_moment_arm_ + Mz1/rotor_moment_constant_ + Tau1);
    controls[0].f2 = 0.25F*(Mx1/quad_lateral_moment_arm_ - My1/quad_longitudinal_moment_arm_ + Mz1/rotor_moment_constant_ + Tau1);
    controls[0].f3 = 0.25F*(Mx1/quad_lateral_moment_arm_ + My1/quad_longitudinal_moment_arm_ - Mz1/rotor_moment_constant_ + Tau1);
    controls[0].f4 = 0.25F*(-Mx1/quad_lateral_moment_arm_ - My1/quad_longitudinal_moment_arm_ - Mz1/rotor_moment_constant_ + Tau1);

    controls[1].f1 = 0.25F*(-Mx2/quad_lateral_moment_arm_ + My2/quad_longitudinal_moment_arm_ + Mz2/rotor_moment_constant_ + Tau2);
    controls[1].f2 = 0.25F*(Mx2/quad_lateral_moment_arm_ - My2/quad_longitudinal_moment_arm_ + Mz2/rotor_moment_constant_ + Tau2);
    controls[1].f3 = 0.25F*(Mx2/quad_lateral_moment_arm_ + My2/quad_longitudinal_moment_arm_ - Mz2/rotor_moment_constant_ + Tau2);
    controls[1].f4 = 0.25F*(-Mx2/quad_lateral_moment_arm_ - My2/quad_longitudinal_moment_arm_ - Mz2/rotor_moment_constant_ + Tau2);
}
