#include "asif_controller.hpp"

namespace asif {
ASIF_Controller::~ASIF_Controller() {
    delete asif_p_;
}

ASIF_Controller::ASIF_Controller(ros::NodeHandle& node_handle) : node_handle_(node_handle), workspace_(workspace), dt_control_(kDefaultControlRate),
    dt_fwd_traj_(kDefaultFwdTrajTimestep),
    quad1_name_(kDefaultQuad1Name), quad2_name_(kDefaultQuad2Name), actuator_topic_(kDefaultQuadActuatorPubTopic),
    gravity_(kDefaultGravity), mass_(kDefaultMass), moment_inertia_x_(kDefaultMomentInertiaX),
    asif_alpha_(kDefaultAsifAlpha), safe_distance_y_(kDefaultSafeDistanceY), safe_distance_z_(kDefaultSafeDistanceZ),
    quad_lateral_moment_arm_(kDefaultQuadLateralMomentArm), quad_longitudinal_moment_arm_(kDefaultQuadLongitudinalMomentArm),
    quad_max_motor_thrust_(kDefaultRotorMaxThrust), rotor_thrust_constant_(kDefaultRotorThrustConstant),
    rotor_moment_constant_(kDefaultRotorMomentConstant), rotor_max_rot_velocity_(kDefaultRotorMaxRotationalVelocity),
    spring_saturation_(kDefaultAsifSpringSaturation), spring_constant_(kDefaultAsifSpringConstant), k_angle_error_(kDefaultAsifAngleErrorGain),
    k_roll_rate_damp_(kDefaultAsifRollRateDamp), k_velocity_damp_(kDefaultAsifVelocityDamp), geom_spring_saturation_(kDefaultGeomSpringSaturation),
    geom_spring_constant_(kDefaultGeomSpringConstant), geom_k_angle_error_(kDefaultGeomAngleErrorGain), geom_k_roll_rate_damp_(kDefaultGeomRollRateDamp),
    geom_k_velocity_damp_(kDefaultGeomVelocityDamp){

    quad1_control_pub_ = node_handle_.advertise<mavros_msgs::ActuatorControl>("/" + quad1_name_ + actuator_topic_, 1);
    quad2_control_pub_ = node_handle_.advertise<mavros_msgs::ActuatorControl>("/" + quad2_name_ + actuator_topic_, 1);
    quad_asif_pub_ = node_handle_.advertise<mav_asif_msgs::ASIF>("/asif/output",1);
    quad_asif_force_moment_pub_ = node_handle_.advertise<mav_asif_msgs::control>("/asif/force_moments",1);
    quad_geom_force_moment_pub_ = node_handle_.advertise<mav_asif_msgs::control>("/desired_controller/force_moments",1);

//    quad1_state_sub_ = node_handle_.subscribe("/" + quad1_name_ + "/mavros/local_position/odom",1, &ASIF_Controller::quad1_state_callback,this);
//    quad2_state_sub_ = node_handle_.subscribe("/" + quad2_name_ + "/mavros/local_position/odom",1, &ASIF_Controller::quad2_state_callback,this);
    quad_state_sub_ = node_handle_.subscribe("/gazebo/model_states",1,&ASIF_Controller::quad_state_callback,this);
    asif_p_ = new asif::ASIF(&workspace_, gravity_, mass_, moment_inertia_x_, spring_saturation_, spring_constant_,
                             asif_alpha_, k_angle_error_, k_roll_rate_damp_, k_velocity_damp_, safe_distance_y_,
                             safe_distance_z_, quad_max_motor_thrust_, quad_lateral_moment_arm_, quad_longitudinal_moment_arm_,
                             rotor_moment_constant_);
}

void ASIF_Controller::quad_state_callback(const gazebo_msgs::ModelStates::ConstPtr& msg)
{
    if (!sensor_data_initialized_) {
        sensor_data_initialized_ = true;
        if (!(msg->name[2].compare("iris0"))){
            iris0_idx_ = 2;
            iris1_idx_ = 3;
        } else {
            iris0_idx_ = 3;
            iris1_idx_ = 2;
        }
    }

    // ----------------------------- iris0_idx_ --------------------------
    Eigen::Vector3d quad1_pos_enu(msg->pose[iris0_idx_].position.x,
                                  msg->pose[iris0_idx_].position.y,
                                  msg->pose[iris0_idx_].position.z);


    Eigen::Vector3d quad1_pos_ned(mavros::ftf::transform_frame_enu_ned<Eigen::Vector3d>(quad1_pos_enu));


    quads_[0].x = static_cast<float>(quad1_pos_ned.x());
    quads_[0].y = static_cast<float>(quad1_pos_ned.y());
    quads_[0].z = static_cast<float>(quad1_pos_ned.z());

    Eigen::Quaterniond q1_enu_baselink((msg->pose[iris0_idx_].orientation.w),
                                       (msg->pose[iris0_idx_].orientation.x),
                                       (msg->pose[iris0_idx_].orientation.y),
                                       (msg->pose[iris0_idx_].orientation.z));

    Eigen::Quaterniond q1_ned_baselink(mavros::ftf::transform_orientation_enu_ned<Eigen::Quaterniond>(q1_enu_baselink));
    Eigen::Quaterniond q1_ned_aircraft(mavros::ftf::transform_orientation_baselink_aircraft<Eigen::Quaterniond>(q1_ned_baselink));

    Eigen::Vector3d orient1_euler = quaternion_to_rpy_wrap(q1_ned_aircraft);

    quads_[0].phi = static_cast<float>(orient1_euler.x());
    quads_[0].theta = static_cast<float>(orient1_euler.y());
    quads_[0].psi = static_cast<float>(orient1_euler.z());

    Eigen::Vector3d quad1_vel_enu(msg->twist[iris0_idx_].linear.x,
                                  msg->twist[iris0_idx_].linear.y,
                                  msg->twist[iris0_idx_].linear.z);

    Eigen::Vector3d quad1_vel_ned(mavros::ftf::transform_frame_enu_ned<Eigen::Vector3d>(quad1_vel_enu));

    quads_[0].vx = static_cast<float>(quad1_vel_ned.x());
    quads_[0].vy = static_cast<float>(quad1_vel_ned.y());
    quads_[0].vz = static_cast<float>(quad1_vel_ned.z());


    Eigen::Vector3d quad1_euler_rates_enu(msg->twist[iris0_idx_].angular.x,
                                          msg->twist[iris0_idx_].angular.y,
                                          msg->twist[iris0_idx_].angular.z);

    Eigen::Vector3d quad1_body_rates_baselink = euler_rates_to_body_rates(quaternion_to_rpy_wrap(q1_enu_baselink),quad1_euler_rates_enu);
    Eigen::Vector3d quad1_body_rates_aircraft(mavros::ftf::transform_frame_baselink_aircraft<Eigen::Vector3d>(quad1_body_rates_baselink));

//    Eigen::Vector3d quad1_body_rates_baselink(mavros::ftf::transform_frame_enu_baselink<Eigen::Vector3d>(quad1_euler_rates_enu,q1_enu_baselink));
//    Eigen::Vector3d quad1_body_rates_aircraft(mavros::ftf::transform_frame_baselink_aircraft<Eigen::Vector3d>(quad1_body_rates_baselink));

    quads_[0].wx = static_cast<float>(-quad1_body_rates_aircraft.y());
    quads_[0].wy = static_cast<float>(quad1_body_rates_aircraft.x());
    quads_[0].wz = static_cast<float>(quad1_body_rates_aircraft.z());


    // ----------------------------- iris1_idx_ --------------------------
    Eigen::Vector3d quad2_pos_enu(msg->pose[iris1_idx_].position.x,
                                  msg->pose[iris1_idx_].position.y,
                                  msg->pose[iris1_idx_].position.z);


    Eigen::Vector3d quad2_pos_ned(mavros::ftf::transform_frame_enu_ned<Eigen::Vector3d>(quad2_pos_enu));


    quads_[1].x = static_cast<float>(quad2_pos_ned.x());
    quads_[1].y = static_cast<float>(quad2_pos_ned.y());
    quads_[1].z = static_cast<float>(quad2_pos_ned.z());

    Eigen::Quaterniond q2_enu_baselink((msg->pose[iris1_idx_].orientation.w),
                                       (msg->pose[iris1_idx_].orientation.x),
                                       (msg->pose[iris1_idx_].orientation.y),
                                       (msg->pose[iris1_idx_].orientation.z));

    Eigen::Quaterniond q2_ned_baselink(mavros::ftf::transform_orientation_enu_ned<Eigen::Quaterniond>(q2_enu_baselink));
    Eigen::Quaterniond q2_ned_aircraft(mavros::ftf::transform_orientation_baselink_aircraft<Eigen::Quaterniond>(q2_ned_baselink));

    Eigen::Vector3d orient2_euler = quaternion_to_rpy_wrap(q2_ned_aircraft);

    quads_[1].phi = static_cast<float>(orient2_euler.x());
    quads_[1].theta = static_cast<float>(orient2_euler.y());
    quads_[1].psi = static_cast<float>(orient2_euler.z());

    Eigen::Vector3d quad2_vel_enu(msg->twist[iris1_idx_].linear.x,
                                  msg->twist[iris1_idx_].linear.y,
                                  msg->twist[iris1_idx_].linear.z);

    Eigen::Vector3d quad2_vel_ned(mavros::ftf::transform_frame_enu_ned<Eigen::Vector3d>(quad2_vel_enu));

    quads_[1].vx = static_cast<float>(quad2_vel_ned.x());
    quads_[1].vy = static_cast<float>(quad2_vel_ned.y());
    quads_[1].vz = static_cast<float>(quad2_vel_ned.z());


    Eigen::Vector3d quad2_euler_rates_enu(msg->twist[iris1_idx_].angular.x,
                                          msg->twist[iris1_idx_].angular.y,
                                          msg->twist[iris1_idx_].angular.z);

    Eigen::Vector3d quad2_body_rates_baselink = euler_rates_to_body_rates(quaternion_to_rpy_wrap(q2_enu_baselink),quad2_euler_rates_enu);
    Eigen::Vector3d quad2_body_rates_aircraft(mavros::ftf::transform_frame_baselink_aircraft<Eigen::Vector3d>(quad2_body_rates_baselink));

//    Eigen::Vector3d quad2_body_rates_baselink(mavros::ftf::transform_frame_enu_baselink<Eigen::Vector3d>(quad2_euler_rates_enu,q2_enu_baselink));
//    Eigen::Vector3d quad2_body_rates_aircraft(mavros::ftf::transform_frame_baselink_aircraft<Eigen::Vector3d>(quad2_body_rates_baselink));

    quads_[1].wx = static_cast<float>(-quad2_body_rates_aircraft.y());
    quads_[1].wy = static_cast<float>(quad2_body_rates_aircraft.x());
    quads_[1].wz = static_cast<float>(quad2_body_rates_aircraft.z());
}
//
//void ASIF_Controller::quad1_state_callback(const nav_msgs::Odometry::ConstPtr& msg)
//{
//    if (!sensor_data_initialized_) {
//        sensor_data_initialized_ = true;
//    }
//
//    Eigen::Vector3d quad1_pos_enu(msg->pose.pose.position.x,
//                                  msg->pose.pose.position.y,
//                                  msg->pose.pose.position.z);
//
//
//    Eigen::Vector3d quad1_pos_ned(mavros::ftf::transform_frame_enu_ned<Eigen::Vector3d>(quad1_pos_enu));
//
//
//    quads_test_[0].x = static_cast<float>(quad1_pos_ned.x());
//    quads_test_[0].y = static_cast<float>(quad1_pos_ned.y());
//    quads_test_[0].z = static_cast<float>(quad1_pos_ned.z());
//
//    Eigen::Quaterniond q1_enu_baselink((msg->pose.pose.orientation.w),
//                                       (msg->pose.pose.orientation.x),
//                                       (msg->pose.pose.orientation.y),
//                                       (msg->pose.pose.orientation.z));
//
//    Eigen::Quaterniond q1_ned_baselink(mavros::ftf::transform_orientation_enu_ned<Eigen::Quaterniond>(q1_enu_baselink));
//    Eigen::Quaterniond q1_ned_aircraft(mavros::ftf::transform_orientation_baselink_aircraft<Eigen::Quaterniond>(q1_ned_baselink));
//
//    Eigen::Vector3d orient1_euler = quaternion_to_rpy_wrap(q1_ned_aircraft);
//
//    quads_test_[0].phi = static_cast<float>(orient1_euler.x());
//    quads_test_[0].theta = static_cast<float>(orient1_euler.y());
//    quads_test_[0].psi = static_cast<float>(orient1_euler.z());
//
//    Eigen::Vector3d quad1_vel_baselink(msg->twist.twist.linear.x,
//                                       msg->twist.twist.linear.y,
//                                       msg->twist.twist.linear.z);
//
//    Eigen::Vector3d quad1_vel_enu(mavros::ftf::transform_frame_baselink_enu<Eigen::Vector3d>(quad1_vel_baselink,q1_enu_baselink));
//    Eigen::Vector3d quad1_vel_ned(mavros::ftf::transform_frame_enu_ned<Eigen::Vector3d>(quad1_vel_enu));
//
//    quads_test_[0].vx = static_cast<float>(quad1_vel_ned.x());
//    quads_test_[0].vy = static_cast<float>(quad1_vel_ned.y());
//    quads_test_[0].vz = static_cast<float>(quad1_vel_ned.z());
//
//
//    Eigen::Vector3d quad1_body_rates_baselink(msg->twist.twist.angular.x,
//                                              msg->twist.twist.angular.y,
//                                              msg->twist.twist.angular.z);
//
//    Eigen::Vector3d quad1_body_rates_aircraft(mavros::ftf::transform_frame_baselink_aircraft<Eigen::Vector3d>(quad1_body_rates_baselink));
//
//    quads_test_[0].wx = static_cast<float>(quad1_body_rates_aircraft.x());
//    quads_test_[0].wy = static_cast<float>(quad1_body_rates_aircraft.y());
//    quads_test_[0].wz = static_cast<float>(quad1_body_rates_aircraft.z());
//}
//
//void ASIF_Controller::quad2_state_callback(const nav_msgs::Odometry::ConstPtr &msg)
//{
//
//    Eigen::Vector3d quad2_pos_enu(msg->pose.pose.position.x + 1.0,
//                                  msg->pose.pose.position.y,
//                                  msg->pose.pose.position.z);
//
//    Eigen::Vector3d quad2_pos_ned(mavros::ftf::transform_frame_enu_ned<Eigen::Vector3d>(quad2_pos_enu));
//
//
//    quads_test_[1].x = static_cast<float>(quad2_pos_ned.x());
//    quads_test_[1].y = static_cast<float>(quad2_pos_ned.y());
//    quads_test_[1].z = static_cast<float>(quad2_pos_ned.z());
//
//
//
//    Eigen::Quaterniond q2_enu_baselink((msg->pose.pose.orientation.w),
//                                       (msg->pose.pose.orientation.x),
//                                       (msg->pose.pose.orientation.y),
//                                       (msg->pose.pose.orientation.z));
//
//    Eigen::Quaterniond q2_ned_baselink(mavros::ftf::transform_orientation_enu_ned<Eigen::Quaterniond>(q2_enu_baselink));
//    Eigen::Quaterniond q2_ned_aircraft(mavros::ftf::transform_orientation_baselink_aircraft<Eigen::Quaterniond>(q2_ned_baselink));
//
//    Eigen::Vector3d orient2_euler = quaternion_to_rpy_wrap(q2_ned_aircraft);
//
//    quads_test_[1].phi = static_cast<float>(orient2_euler.x());
//    quads_test_[1].theta = static_cast<float>(orient2_euler.y());
//    quads_test_[1].psi = static_cast<float>(orient2_euler.z());
//
//    Eigen::Vector3d quad2_vel_baselink(msg->twist.twist.linear.x,
//                                       msg->twist.twist.linear.y,
//                                       msg->twist.twist.linear.z);
//    Eigen::Vector3d quad2_vel_enu(mavros::ftf::transform_frame_baselink_enu<Eigen::Vector3d>(quad2_vel_baselink,q2_enu_baselink));
//    Eigen::Vector3d quad2_vel_ned(mavros::ftf::transform_frame_enu_ned<Eigen::Vector3d>(quad2_vel_enu));
//
//    quads_test_[1].vx = static_cast<float>(quad2_vel_ned.x());
//    quads_test_[1].vy = static_cast<float>(quad2_vel_ned.y());
//    quads_test_[1].vz = static_cast<float>(quad2_vel_ned.z());
//
//    Eigen::Vector3d quad2_body_rates_baselink(msg->twist.twist.angular.x,
//                                              msg->twist.twist.angular.y,
//                                              msg->twist.twist.angular.z);
//
//    Eigen::Vector3d quad2_body_rates_aircraft(mavros::ftf::transform_frame_baselink_aircraft<Eigen::Vector3d>(quad2_body_rates_baselink));
//
//    quads_test_[1].wx = static_cast<float>(quad2_body_rates_aircraft.x());
//    quads_test_[1].wy = static_cast<float>(quad2_body_rates_aircraft.y());
//    quads_test_[1].wz = static_cast<float>(quad2_body_rates_aircraft.z());
//}

void
ASIF_Controller::geometric_controller(const asif::quad_states_s pos_des[2], asif::control_s controls[2]) {
    float f_cmd_x1 = geom_spring_constant_ * tanhf(geom_spring_saturation_ * (pos_des[0].x - quads_[0].x)) -
            geom_k_velocity_damp_ * quads_[0].vx;
    float f_cmd_x2 = geom_spring_constant_ * tanhf(geom_spring_saturation_ * (pos_des[1].x - quads_[1].x)) -
            geom_k_velocity_damp_ * quads_[1].vx;
    float f_cmd_y1 = geom_spring_constant_ * tanhf(geom_spring_saturation_ * (pos_des[0].y - quads_[0].y)) -
            geom_k_velocity_damp_ * quads_[0].vy;
    float f_cmd_y2 = geom_spring_constant_ * tanhf(geom_spring_saturation_ * (pos_des[1].y - quads_[1].y)) -
            geom_k_velocity_damp_ * quads_[1].vy;
    float f_cmd_z1 = geom_spring_constant_ * tanhf(geom_spring_saturation_ * (pos_des[0].z - quads_[0].z)) - mass_ * gravity_ -
            geom_k_velocity_damp_ * quads_[0].vz;
    float f_cmd_z2 = geom_spring_constant_ * tanhf(geom_spring_saturation_ * (pos_des[1].z - quads_[1].z)) - mass_ * gravity_ -
            geom_k_velocity_damp_ * quads_[1].vz;

    Eigen::Vector3f thrust_dir1 {f_cmd_x1, f_cmd_y1, f_cmd_z1};
    Eigen::Vector3f thrust_dir2 {f_cmd_x2, f_cmd_y2, f_cmd_z2};
    thrust_dir1.normalize();
    thrust_dir2.normalize();
    float psi1_cmd = 0;
    float psi2_cmd = 0;

    Eigen::Vector3d iterm_thrust_frame1_1 = mavros::ftf::transform_frame_ned_aircraft
            <Eigen::Vector3d>(thrust_dir1.cast<double>(),mavros::ftf::quaternion_from_rpy(0,0,(double)psi1_cmd));
    Eigen::Vector3d iterm_thrust_frame1_2 = mavros::ftf::transform_frame_ned_aircraft
            <Eigen::Vector3d>(thrust_dir2.cast<double>(), mavros::ftf::quaternion_from_rpy(0,0,(double)psi2_cmd));


    float theta1_cmd = atan2f(-iterm_thrust_frame1_1.cast<float>().x(),-iterm_thrust_frame1_1.cast<float>().z());
    float theta2_cmd = atan2f(-iterm_thrust_frame1_2.cast<float>().x(),-iterm_thrust_frame1_2.cast<float>().z());

    Eigen::Vector3d iterm_thrust_frame2_1 = mavros::ftf::transform_frame_ned_aircraft
            <Eigen::Vector3d>(iterm_thrust_frame1_1,mavros::ftf::quaternion_from_rpy(0,(double)theta1_cmd,0));
    Eigen::Vector3d iterm_thrust_frame2_2 = mavros::ftf::transform_frame_ned_aircraft
            <Eigen::Vector3d>(iterm_thrust_frame1_2, mavros::ftf::quaternion_from_rpy(0,(double)theta2_cmd,0));


    float phi1_cmd = atan2f(iterm_thrust_frame2_1.cast<float>().y(),-iterm_thrust_frame2_1.cast<float>().z());
    float phi2_cmd = atan2f(iterm_thrust_frame2_2.cast<float>().y(),-iterm_thrust_frame2_2.cast<float>().z());

    float Tau1 = - f_cmd_x1 * cosf(quads_[0].phi) * sinf(quads_[0].theta) + f_cmd_y1 * sinf(quads_[0].phi)
            - f_cmd_z1 * cosf(quads_[0].phi) * cosf(quads_[0].theta);
    float Tau2 = - f_cmd_x2 * cosf(quads_[1].phi) * sinf(quads_[1].theta) + f_cmd_y2 * sinf(quads_[1].phi)
            - f_cmd_z2 * cosf(quads_[1].phi) * cosf(quads_[1].theta);

    float Mx1 = -geom_k_angle_error_ * (quads_[0].phi - phi1_cmd) - geom_k_roll_rate_damp_ * quads_[0].wx;
    float Mx2 = -geom_k_angle_error_ * (quads_[1].phi - phi2_cmd) - geom_k_roll_rate_damp_ * quads_[1].wx;

    float My1 = -geom_k_angle_error_ * (quads_[0].theta - theta1_cmd) - geom_k_roll_rate_damp_ * quads_[0].wy;
    float My2 = -geom_k_angle_error_ * (quads_[1].theta - theta2_cmd) - geom_k_roll_rate_damp_ * quads_[1].wy;

    float Mz1 = -geom_k_angle_error_ * (quads_[0].psi - psi1_cmd) - geom_k_roll_rate_damp_ * quads_[0].wz;
    float Mz2 = -geom_k_angle_error_ * (quads_[1].psi - psi2_cmd) - geom_k_roll_rate_damp_ * quads_[1].wz;

    quad_geom_force_moment_msg_.moment1 = Mx1;
    quad_geom_force_moment_msg_.moment2 = Mx2;
    quad_geom_force_moment_msg_.thrust1 = Tau1;
    quad_geom_force_moment_msg_.thrust2 = Tau2;

    quad_geom_force_moment_pub_.publish(quad_geom_force_moment_msg_);

//    printf("\n\n--------------------gazebo_model_states---------------------\n");
//    printf("x1: %f, y1: %f, z1: %f\n",quads_[0].x,quads_[0].y,quads_[0].z);
//    printf("x2: %f, y2: %f, z2: %f\n",quads_[1].x,quads_[1].y,quads_[1].z);
//    printf("vx1: %f, vy1: %f, z1: %f\n",quads_[0].vx,quads_[0].vy,quads_[0].vz);
//    printf("vx2: %f, vy2: %f, z2: %f\n",quads_[1].vx,quads_[1].vy,quads_[1].vz);
//    printf("phi1: %f, theta1: %f, psi1: %f\n", quads_[0].phi, quads_[0].theta, quads_[0].psi);
//    printf("phi2: %f, theta2: %f, psi2: %f\n", quads_[1].phi, quads_[1].theta, quads_[1].psi);
//    printf("wx1: %f, wy1: %f, wz1: %f\n", quads_[0].wx, quads_[0].wy, quads_[0].wz);
//    printf("wx2: %f, wy2: %f, wz2: %f\n", quads_[1].wx, quads_[1].wy, quads_[1].wz);
//    printf("--------------------local_odom_states---------------------\n");
//    printf("x1: %f, y1: %f, z1: %f\n",quads_test_[0].x,quads_test_[0].y,quads_test_[0].z);
//    printf("x2: %f, y2: %f, z2: %f\n",quads_test_[1].x,quads_test_[1].y,quads_test_[1].z);
//    printf("vx1: %f, vy1: %f, z1: %f\n",quads_test_[0].vx,quads_test_[0].vy,quads_test_[0].vz);
//    printf("vx2: %f, vy2: %f, z2: %f\n",quads_test_[1].vx,quads_test_[1].vy,quads_test_[1].vz);
//    printf("phi1: %f, theta1: %f, psi1: %f\n", quads_test_[0].phi, quads_test_[0].theta, quads_test_[0].psi);
//    printf("phi2: %f, theta2: %f, psi2: %f\n", quads_test_[1].phi, quads_test_[1].theta, quads_test_[1].psi);
//    printf("wx1: %f, wy1: %f, wz1: %f\n", quads_test_[0].wx, quads_test_[0].wy, quads_test_[0].wz);
//    printf("wx2: %f, wy2: %f, wz2: %f\n", quads_test_[1].wx, quads_test_[1].wy, quads_test_[1].wz);

    controls[0].f1 = 0.25F*(-Mx1/quad_lateral_moment_arm_ + My1/quad_longitudinal_moment_arm_ + Mz1/rotor_moment_constant_ + Tau1);
    controls[0].f2 = 0.25F*(Mx1/quad_lateral_moment_arm_ - My1/quad_longitudinal_moment_arm_ + Mz1/rotor_moment_constant_ + Tau1);
    controls[0].f3 = 0.25F*(Mx1/quad_lateral_moment_arm_ + My1/quad_longitudinal_moment_arm_ - Mz1/rotor_moment_constant_ + Tau1);
    controls[0].f4 = 0.25F*(-Mx1/quad_lateral_moment_arm_ - My1/quad_longitudinal_moment_arm_ - Mz1/rotor_moment_constant_ + Tau1);

    controls[1].f1 = 0.25F*(-Mx2/quad_lateral_moment_arm_ + My2/quad_longitudinal_moment_arm_ + Mz2/rotor_moment_constant_ + Tau2);
    controls[1].f2 = 0.25F*(Mx2/quad_lateral_moment_arm_ - My2/quad_longitudinal_moment_arm_ + Mz2/rotor_moment_constant_ + Tau2);
    controls[1].f3 = 0.25F*(Mx2/quad_lateral_moment_arm_ + My2/quad_longitudinal_moment_arm_ - Mz2/rotor_moment_constant_ + Tau2);
    controls[1].f4 = 0.25F*(-Mx2/quad_lateral_moment_arm_ - My2/quad_longitudinal_moment_arm_ - Mz2/rotor_moment_constant_ + Tau2);
}

int ASIF_Controller::asif_filter(asif::control_s controls[2])
{
    matrix::Vector<float,20> phi_max;
    matrix::Vector<float,4> force_moments;
    float psi;
    int phi_max_time;
    if (asif_p_->QP(quads_, controls,phi_max, psi, phi_max_time,force_moments)) {
        quad_asif_msg_.psi = psi;
        quad_asif_msg_.maxTime = phi_max_time;
        for (int i = 0; i < 20; i ++){
            quad_asif_msg_.phiMax[i] = phi_max(i);
        }
        quad_asif_force_moment_msg_.thrust1 = force_moments(0);
        quad_asif_force_moment_msg_.thrust2 = force_moments(1);
        quad_asif_force_moment_msg_.moment1 = force_moments(2);
        quad_asif_force_moment_msg_.moment2 = force_moments(3);

        quad_asif_pub_.publish(quad_asif_msg_);
        quad_asif_force_moment_pub_.publish(quad_asif_force_moment_msg_);
        return 1;
    }
    return 0;
}

void ASIF_Controller::backup_control(asif::control_s controls[2])
{
    asif_p_->backupControl(quads_,controls);
}

float ASIF_Controller::motor_force_to_mixer_input(float force) const
{
    if (force < 0) {
        force = 0;
    } else if (force > quad_max_motor_thrust_) {
        force = quad_max_motor_thrust_;
    }
    float rot_velocity = sqrtf(force / rotor_thrust_constant_)*1.000661664F;
    float mixer_input = (rot_velocity - 100.0F) / 1000.0F;
    return mixer_input;
}

void ASIF_Controller::update_actuator_quads(asif::control_s controls[2])
{
    quad1_actuator_msg_.group_mix = 0;
    quad2_actuator_msg_.group_mix = 0;

    quad1_actuator_msg_.controls[0] = motor_force_to_mixer_input(controls[0].f1);
    quad1_actuator_msg_.controls[1] = motor_force_to_mixer_input(controls[0].f2);
    quad1_actuator_msg_.controls[2] = motor_force_to_mixer_input(controls[0].f3);
    quad1_actuator_msg_.controls[3] = motor_force_to_mixer_input(controls[0].f4);
    quad2_actuator_msg_.controls[0] = motor_force_to_mixer_input(controls[1].f1);
    quad2_actuator_msg_.controls[1] = motor_force_to_mixer_input(controls[1].f2);
    quad2_actuator_msg_.controls[2] = motor_force_to_mixer_input(controls[1].f3);
    quad2_actuator_msg_.controls[3] = motor_force_to_mixer_input(controls[1].f4);
}

void ASIF_Controller::actuatorTimerCallback(const ros::TimerEvent &)
{
    mavros_msgs::ActuatorControl quad1_actuator_msg_tmp = quad1_actuator_msg_;
    mavros_msgs::ActuatorControl quad2_actuator_msg_tmp = quad2_actuator_msg_;

    quad1_control_pub_.publish(quad1_actuator_msg_);
    quad2_control_pub_.publish(quad2_actuator_msg_);

    quad1_actuator_msg_ = quad1_actuator_msg_tmp;
    quad2_actuator_msg_ = quad2_actuator_msg_tmp;
}

void ASIF_Controller::getQuadsMsg(mavros_msgs::ActuatorControl &quad1_actuator_msg, mavros_msgs::ActuatorControl &quad2_actuator_msg)
{
    quad1_actuator_msg = quad1_actuator_msg_;
    quad2_actuator_msg = quad2_actuator_msg_;
}

Eigen::Vector3d ASIF_Controller::quaternion_to_rpy_wrap(const Eigen::Quaterniond& q)
{
    Eigen::Vector3d rpy;
    double roll = atan2(2 * (q.w() * q.x() + q.y() * q.z()),1 - 2 * (pow(q.x(),2) + pow(q.y(),2)));
    double pitch = asin(2 * (q.w() * q.y() - q.z() * q.x()));
    double yaw = atan2(2 * (q.w() * q.z() + q.x() * q.y()),1 - 2 *(pow(q.y(),2) + pow(q.z(),2)));

    rpy << roll,
            pitch,
            yaw;

    return rpy;
}

Eigen::Vector3d ASIF_Controller::euler_rates_to_body_rates(const Eigen::Vector3d orientation,const Eigen::Vector3d& euler_rates)
{
    Eigen::Vector3d body_rates;

    body_rates << euler_rates.x() - sin(orientation.y())*euler_rates.z(),
                  cos(orientation.x())*euler_rates.y() + sin(orientation.x())*cos(orientation.y())*euler_rates.z(),
                  -sin(orientation.x())*euler_rates.y() + cos(orientation.x())*cos(orientation.y())*euler_rates.z();

    return body_rates;
}

bool ASIF_Controller::sensors_initialized()
{
    return sensor_data_initialized_;
}


}


