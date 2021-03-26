#pragma once

#include "osqp.h"
#include "matrix/matrix/math.hpp"

namespace asif {

typedef struct {
    float x;
    float y;
    float z;
    float vx;
    float vy;
    float vz;
    float phi;
    float theta;
    float psi;
    float wx;
    float wy;
    float wz;
} quad_states_s;

typedef struct {
    float f1;
    float f2;
    float f3;
    float f4;
} control_s;

typedef float states_t[10];
typedef matrix::Vector<float,20> embed_state_t;
typedef float disturbance_corners_t[6][64];
typedef float state_space_corners_t[10][1024];
typedef float alpha_beta_corners_t[10][4];

#define DT_BACKUP 0.01
#define T_BACKUP 5
#define BARRIER_CONSTANT_BACKUP 0.02
#define SAFE_RADIUS 0.3
#define NUM_BACKUP_POINTS 5
#define BACKUP_STEPS static_cast<int>(T_BACKUP/DT_BACKUP + 1)

class ASIF {
public:
    ASIF(OSQPWorkspace *workspace, float gravity, float mass, float moment_inertia_x, float spring_saturation,
         float spring_constant, float asif_alpha, float k_angle_error, float k_roll_rate_damp, float k_velocity_damp,
         float safe_distance_y, float safe_distance_z, float quad_max_thrust, float quad_lateral_moment_arm,
         float quad_longitudinal_moment_arm, float rotor_moment_constant);
    ~ASIF();

    int QP(const quad_states_s *quads, control_s *controls, matrix::Vector<float,20> &phi_max, float &psi, int &phi_max_time, matrix::Vector<float,4> &force_moments);

//    void dynamics(const quad_states_s states[2], const control_s controls[2], quad_states_s states_dot[2]);
    void backupControl(const quad_states_s states[2], control_s controls[2]);

private:
    float gravity_{};
    float mass_{};
    float moment_inertia_x_{};
    float spring_saturation_{};
    float spring_constant_{};
    float asif_alpha_{};
    float k_angle_error_{};
    float k_roll_rate_damp_{};
    float k_velocity_damp_{};
    float safe_distance_y_{};
    float safe_distance_z_{};
    float quad_lateral_moment_arm_{};
    float quad_longitudinal_moment_arm_{};
    float rotor_moment_constant_{};
    float quad_max_thrust_{};

    float w_min_[6] = {-0.00001,-0.00001,-0.0002,-0.0002,-0.000000001,-0.000000001};
    float w_max_[6] = {0.00001,0.00001,0.0002,0.0002,0.000000001,0.000000001};
    disturbance_corners_t w_corners_{};

    float xbar_[10]{};

    float Pf_[10][10]{};
    matrix::Matrix<float,10,10> matrix_Pf_{};
    matrix::Matrix<float,10,10> matrix_Sf_{};

    static constexpr float barrier_constant_backup_{BARRIER_CONSTANT_BACKUP};
    static constexpr float safe_radius_{SAFE_RADIUS};

    static constexpr float dt_backup_{DT_BACKUP};
    static constexpr float T_backup_{T_BACKUP};
    static constexpr int backup_barrier_steps_{NUM_BACKUP_POINTS};
    static constexpr int backup_steps_{BACKUP_STEPS};
    int max_time_ind_backup_{-1};
    int min_time_ind_unsafe_{-1};

    bool QP_initialized{false};
    matrix::Matrix<float,10,20> derivative_corners_[1024];
    matrix::Matrix<float,10,20> derivative_corners_unsafe_[4];

    OSQPWorkspace *osqp_workspace_;
    embed_state_t *phi_;

    matrix::Vector<float,10> decomp(const states_t x, const states_t xh, const float w[6]) const;
    embed_state_t embeddingDynamics(embed_state_t phi);
    matrix::Matrix<float,10,10> jacobian(states_t x);
    void phiE(states_t x0, states_t xh0);
    void hsamCorners(state_space_corners_t corners, float *hs);
    void hsamCornersUnsafe(alpha_beta_corners_t corners, float *hs);
    static void getCorners(const float *w_min, const float *w_max, disturbance_corners_t corners);
    static void getCornersUnsafe(states_t x_min, states_t x_max, alpha_beta_corners_t corners);
    static void getCorners(states_t x_min, states_t x_max, state_space_corners_t corners);
    void getDerivativeCorners();
    void backupTrajectoryDistribution(int traj_t, int backup_steps, int *out);
    void hSamSoftMin(embed_state_t phi, float *hs_min);
    void hSamSoftMinUnsafe(embed_state_t phi, float *hs_min);
    matrix::Matrix<float,1,20> hSamSoftMinGrad(embed_state_t phi);
    matrix::Matrix<float,1,20> hSamSoftMinGradUnsafe(embed_state_t phi);
    
};
}
