#pragma once

#include "osqp.h"
#include "matrix/matrix/math.hpp"

namespace asif {

typedef struct {
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double phi;
    double theta;
    double psi;
    double wx;
    double wy;
    double wz;
} quad_states_s;

typedef struct {
    double f1;
    double f2;
    double f3;
    double f4;
} control_s;

typedef double states_t[10];
typedef matrix::Vector<double,20> embed_state_t;
typedef double disturbance_corners_t[6][64];
typedef double state_space_corners_t[10][1024];
typedef double alpha_beta_corners_t[10][4];

class ASIF {
public:
    ASIF(OSQPWorkspace *workspace, double gravity, double mass, double moment_inertia_x, double spring_saturation,
         double spring_constant, double asif_alpha, double k_angle_error, double k_roll_rate_damp, double k_velocity_damp,
         double safe_distance_y, double safe_distance_z, double quad_max_thrust, double quad_lateral_moment_arm,
         double quad_longitudinal_moment_arm, double rotor_moment_constant);
    ~ASIF();

    int QP(const quad_states_s *quads, control_s *controls);

//    void dynamics(const quad_states_s states[2], const control_s controls[2], quad_states_s states_dot[2]);
    void backupControl(const quad_states_s states[2], control_s controls[2]);

private:
    double gravity_;
    double mass_;
    double moment_inertia_x_;
    double spring_saturation_;
    double spring_constant_;
    double asif_alpha_;
    double k_angle_error_;
    double k_roll_rate_damp_;
    double k_velocity_damp_;
    double safe_distance_y_;
    double safe_distance_z_;
    double quad_lateral_moment_arm_;
    double quad_longitudinal_moment_arm_;
    double rotor_moment_constant_;
    double quad_max_thrust_;

    double w_min_[6] = {-0.000007,-0.000007,-0.000007,-0.000007,-0.000007,-0.000007};
    double w_max_[6] = {0.000007,0.000007,0.000007,0.000007,0.000007,0.000007};
    disturbance_corners_t w_corners_;

    double xbar_[10];

    double Pf_[10][10];
    matrix::Matrix<double,10,10> matrix_Pf_;
    matrix::Matrix<double,10,10> matrix_Sf_;

    double barrier_constant_backup_;
    double barrier_constant_unsafe_;

    double dt_backup_ = 0.01;
    double T_backup_ = 0.5;
    int backup_steps_;
    int backup_barrier_steps_ = 6;
    int max_time_ind_backup_;
    int min_time_ind_unsafe_;

    bool QP_initialized = false;
    matrix::Matrix<double,10,20> derivative_corners_[1024];
    matrix::Matrix<double,10,20> derivative_corners_unsafe_[4];

    OSQPWorkspace *osqp_workspace_;
    embed_state_t *phi_;

    matrix::Vector<double,10> decomp(const states_t x, const states_t xh, const double w[6]) const;
    embed_state_t embeddingDynamics(embed_state_t phi);
    matrix::Matrix<double,10,10> jacobian(states_t x);
    void phiE(states_t x0, states_t xh0);
    void hsamCorners(state_space_corners_t corners, double *hs);
    void hsamCornersUnsafe(alpha_beta_corners_t corners, double *hs);
    static void getCorners(const double *w_min, const double *w_max, disturbance_corners_t corners);
    static void getCornersUnsafe(states_t x_min, states_t x_max, alpha_beta_corners_t corners);
    static void getCorners(states_t x_min, states_t x_max, state_space_corners_t corners);
    void getDerivativeCorners();
    void backupTrajectoryDistribution(int traj_t, int backup_steps, int *out);
    void hSamSoftMin(embed_state_t phi, double *hs_min);
    void hSamSoftMinUnsafe(embed_state_t phi, double *hs_min);
    matrix::Matrix<double,1,20> hSamSoftMinGrad(embed_state_t phi);
    matrix::Matrix<double,1,20> hSamSoftMinGradUnsafe(embed_state_t phi);
    
};
}
