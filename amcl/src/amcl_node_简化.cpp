/*
 *  Copyright (c) 2008, Willow Garage, Inc.
 *  All rights reserved.
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* Author: Brian Gerkey */

#include <algorithm>
#include <vector>
#include <map>
#include <cmath>
#include <memory>

#include <boost/bind.hpp>
#include <boost/thread/mutex.hpp>

// Signal handling
#include <signal.h>

#include "amcl/map/map.h"
#include "amcl/pf/pf.h"
#include "amcl/sensors/amcl_odom.h"
#include "amcl/sensors/amcl_laser.h"
#include "portable_utils.hpp"

#include "ros/assert.h"

// roscpp
#include "ros/ros.h"

// Messages that I need
#include "sensor_msgs/LaserScan.h"
#include "geometry_msgs/PoseWithCovarianceStamped.h"
#include "geometry_msgs/PoseArray.h"
#include "geometry_msgs/Pose.h"
#include "geometry_msgs/PoseStamped.h"
#include "nav_msgs/GetMap.h"
#include "nav_msgs/SetMap.h"
#include "std_srvs/Empty.h"

// For transform support
#include "tf2/LinearMath/Transform.h"
#include "tf2/convert.h"
#include "tf2/utils.h"
#include "tf2_geometry_msgs/tf2_geometry_msgs.h"
#include "tf2_ros/buffer.h"
#include "tf2_ros/message_filter.h"
#include "tf2_ros/transform_broadcaster.h"
#include "tf2_ros/transform_listener.h"
#include "message_filters/subscriber.h"

// Dynamic_reconfigure
#include "dynamic_reconfigure/server.h"
#include "amcl/AMCLConfig.h"

// Allows AMCL to run from bag file
#include <rosbag/bag.h>
#include <rosbag/view.h>
#include <boost/foreach.hpp>

// For monitoring the estimator
#include <diagnostic_updater/diagnostic_updater.h>

#define NEW_UNIFORM_SAMPLING 1

using namespace amcl;

// Pose hypothesis
typedef struct
{
    // Total weight (weights sum to 1)
    double weight;

    // Mean of pose esimate
    pf_vector_t pf_pose_mean;

    // Covariance of pose estimate
    pf_matrix_t pf_pose_cov;

} amcl_hyp_t;

static double
normalize(double z)
{
    return atan2(sin(z), cos(z));
}
static double
angle_diff(double a, double b)
{
    double d1, d2;
    a = normalize(a);
    b = normalize(b);
    d1 = a - b;
    d2 = 2 * M_PI - fabs(d1);
    if (d1 > 0)
        d2 *= -1.0;
    if (fabs(d1) < fabs(d2))
        return (d1);
    else
        return (d2);
}

static const std::string scan_topic_ = "scan";

/* This function is only useful to have the whole code work
 * with old rosbags that have trailing slashes for their frames
 */
inline std::string stripSlash(const std::string &in)
{
    std::string out = in;
    if ((!in.empty()) && (in[0] == '/'))
        out.erase(0, 1);
    return out;
}

class AmclNode
{
public:
    AmclNode();
    ~AmclNode();

    /**
     * @brief Uses TF and LaserScan messages from bag file to drive AMCL instead
     * @param in_bag_fn input bagfile
     * @param trigger_global_localization whether to trigger global localization
     * before starting to process the bagfile
     */
    void runFromBag(const std::string &in_bag_fn, bool trigger_global_localization = false);

    int process();
    void savePoseToServer();

private:
    std::shared_ptr<tf2_ros::TransformBroadcaster> tfb_;
    std::shared_ptr<tf2_ros::TransformListener> tfl_;
    std::shared_ptr<tf2_ros::Buffer> tf_;

    bool sent_first_transform_;

    tf2::Transform latest_tf_;
    bool latest_tf_valid_;

    // Pose-generating function used to uniformly distribute particles over
    // the map
    static pf_vector_t uniformPoseGenerator(void *arg);
#if NEW_UNIFORM_SAMPLING
    static std::vector<std::pair<int, int>> free_space_indices; // class static 成员声明
#endif
    // Callbacks					回调函数声明
    bool globalLocalizationCallback(std_srvs::Empty::Request &req,
                                    std_srvs::Empty::Response &res);
    bool nomotionUpdateCallback(std_srvs::Empty::Request &req,
                                std_srvs::Empty::Response &res);
    bool setMapCallback(nav_msgs::SetMap::Request &req,
                        nav_msgs::SetMap::Response &res);

    void laserReceived(const sensor_msgs::LaserScanConstPtr &laser_scan);
    void initialPoseReceived(const geometry_msgs::PoseWithCovarianceStampedConstPtr &msg);
    void handleInitialPoseMessage(const geometry_msgs::PoseWithCovarianceStamped &msg);
    void mapReceived(const nav_msgs::OccupancyGridConstPtr &msg);

    void handleMapMessage(const nav_msgs::OccupancyGrid &msg);
    void freeMapDependentMemory();
    map_t *convertMap(const nav_msgs::OccupancyGrid &map_msg); // 把订阅的map消息转换为map_t类型的地图
    void updatePoseFromServer();
    void applyInitialPose();

    // parameter for which odom to use
    std::string odom_frame_id_; // odom 里程计坐标系

    // paramater to store latest odom pose
    geometry_msgs::PoseStamped latest_odom_pose_; // latest最新的里程计位姿  latest_odom_pose_

    // parameter for which base to use
    std::string base_frame_id_;   // base_link
    std::string global_frame_id_; // map 坐标系

    bool use_map_topic_;
    bool first_map_only_;

    ros::Duration gui_publish_period;
    ros::Time save_pose_last_time;
    ros::Duration save_pose_period;

    geometry_msgs::PoseWithCovarianceStamped last_published_pose; // 上一次发布的位姿last_published_pose;

    map_t *map_;   // AMCL的地图指针
    char *mapdata; // map.data中的具体栅格数值
    int sx, sy;
    double resolution; // 分辨率
                       // 与激光数据有关的 "订阅", 存储激光的vector，存储更新激光数据的vector
    message_filters::Subscriber<sensor_msgs::LaserScan> *laser_scan_sub_;
    tf2_ros::MessageFilter<sensor_msgs::LaserScan> *laser_scan_filter_;
    ros::Subscriber initial_pose_sub_; // 初始位姿 订阅器
    std::vector<AMCLLaser *> lasers_;
    std::vector<bool> lasers_update_;
    std::map<std::string, int> frame_to_laser_;

    // Particle filter			//  粒子滤波器相关的变量
    pf_t *pf_;
    double pf_err_, pf_z_;
    bool pf_init_;
    pf_vector_t pf_odom_pose_;
    double d_thresh_, a_thresh_;
    int resample_interval_;
    int resample_count_;
    double laser_min_range_;
    double laser_max_range_;

    // Nomotion update control
    bool m_force_update; // used to temporarily let amcl update samples even when no motion occurs...

    AMCLOdom *odom_;
    AMCLLaser *laser_;

    ros::Duration cloud_pub_interval;
    ros::Time last_cloud_pub_time;

    // For slowing play-back when reading directly from a bag file
    ros::WallDuration bag_scan_period_;

    void requestMap();

    // Helper to get odometric pose from transform system
    bool getOdomPose(geometry_msgs::PoseStamped &pose,
                     double &x, double &y, double &yaw,
                     const ros::Time &t, const std::string &f);

    // time for tolerance on the published transform,
    // basically defines how long a map->odom transform is good for
    ros::Duration transform_tolerance_;

    ros::NodeHandle nh_; // 节点，订阅和发布器
    ros::NodeHandle private_nh_;
    ros::Publisher pose_pub_;
    ros::Publisher particlecloud_pub_;
    ros::ServiceServer global_loc_srv_;
    ros::ServiceServer nomotion_update_srv_; // to let amcl update samples without requiring motion
    ros::ServiceServer set_map_srv_;
    ros::Subscriber initial_pose_sub_old_;
    ros::Subscriber map_sub_;

    diagnostic_updater::Updater diagnosic_updater_; //???
    void standardDeviationDiagnostics(diagnostic_updater::DiagnosticStatusWrapper &diagnostic_status);
    double std_warn_level_x_;
    double std_warn_level_y_;
    double std_warn_level_yaw_;

    amcl_hyp_t *initial_pose_hyp_;
    bool first_map_received_;
    bool first_reconfigure_call_;

    boost::recursive_mutex configuration_mutex_;
    dynamic_reconfigure::Server<amcl::AMCLConfig> *dsrv_;
    amcl::AMCLConfig default_config_;
    ros::Timer check_laser_timer_;

    int max_beams_, min_particles_, max_particles_;
    double alpha1_, alpha2_, alpha3_, alpha4_, alpha5_;
    double alpha_slow_, alpha_fast_;
    double z_hit_, z_short_, z_max_, z_rand_, sigma_hit_, lambda_short_;
    // beam skip related params
    bool do_beamskip_;
    double beam_skip_distance_, beam_skip_threshold_, beam_skip_error_threshold_;
    double laser_likelihood_max_dist_;
    odom_model_t odom_model_type_;
    double init_pose_[3];
    double init_cov_[3];
    laser_model_t laser_model_type_;
    bool tf_broadcast_;
    bool force_update_after_initialpose_;
    bool force_update_after_set_map_;
    bool selective_resampling_;

    void reconfigureCB(amcl::AMCLConfig &config, uint32_t level);

    ros::Time last_laser_received_ts_;
    ros::Duration laser_check_interval_;
    void checkLaserReceived(const ros::TimerEvent &event);
};

#if NEW_UNIFORM_SAMPLING
std::vector<std::pair<int, int>> AmclNode::free_space_indices; //  类中的静态数据成员，在外部进行定义和初始化，类内部声明
#endif

#define USAGE "USAGE: amcl"

boost::shared_ptr<AmclNode> amcl_node_ptr; // 主要操作对象

void sigintHandler(int sig)
{
    // Save latest pose as we're shutting down.
    amcl_node_ptr->savePoseToServer();      // amcl节点退出的时候，执行该方法，将某几个数据保存到server(参数服务器中)，此时针对仅仅是amcl节点退出的情况，roscore本身还没有退出，参数服务器还是处于运行状态；
    ros::shutdown();
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "amcl");
    ros::NodeHandle nh;

    // Override default sigint handler
    signal(SIGINT, sigintHandler);

    // Make our node available to sigintHandler
    amcl_node_ptr.reset(new AmclNode());

    if (argc == 1)
    {
        // run using ROS input
        ros::spin();
    }
    else if ((argc >= 3) && (std::string(argv[1]) == "--run-from-bag"))
    {
        if (argc == 3)
        {
            amcl_node_ptr->runFromBag(argv[2]);
        }
        else if ((argc == 4) && (std::string(argv[3]) == "--global-localization"))
        {
            amcl_node_ptr->runFromBag(argv[2], true);
        }
    }

    // Without this, our boost locks are not shut down nicely
    amcl_node_ptr.reset();

    // To quote Morgan, Hooray!
    return (0);
}
// WHOLELINE:   AmclNode()
AmclNode::AmclNode() : // ##step01 构造函数
                       sent_first_transform_(false),
                       latest_tf_valid_(false),
                       map_(NULL), //
                       pf_(NULL),
                       resample_count_(0),
                       odom_(NULL),
                       laser_(NULL),
                       private_nh_("~"),
                       initial_pose_hyp_(NULL),
                       first_map_received_(false),
                       first_reconfigure_call_(true)
{
    boost::recursive_mutex::scoped_lock l(configuration_mutex_);

    // Grab params off the param server
    private_nh_.param("use_map_topic", use_map_topic_, false);
    private_nh_.param("first_map_only", first_map_only_, false);

    double tmp;
    private_nh_.param("gui_publish_rate", tmp, -1.0);
    gui_publish_period = ros::Duration(1.0 / tmp);
    private_nh_.param("save_pose_rate", tmp, 0.5);
    save_pose_period = ros::Duration(1.0 / tmp);

    private_nh_.param("laser_min_range", laser_min_range_, -1.0);
    private_nh_.param("laser_max_range", laser_max_range_, -1.0);
    private_nh_.param("laser_max_beams", max_beams_, 30);
    private_nh_.param("min_particles", min_particles_, 100);
    private_nh_.param("max_particles", max_particles_, 5000);
    private_nh_.param("kld_err", pf_err_, 0.01);    // pf_err_
    private_nh_.param("kld_z", pf_z_, 0.99);        // pf_z_
    private_nh_.param("odom_alpha1", alpha1_, 0.2);
    private_nh_.param("odom_alpha2", alpha2_, 0.2);
    private_nh_.param("odom_alpha3", alpha3_, 0.2);
    private_nh_.param("odom_alpha4", alpha4_, 0.2);
    private_nh_.param("odom_alpha5", alpha5_, 0.2);

    private_nh_.param("do_beamskip", do_beamskip_, false);
    private_nh_.param("beam_skip_distance", beam_skip_distance_, 0.5);
    private_nh_.param("beam_skip_threshold", beam_skip_threshold_, 0.3);
    if (private_nh_.hasParam("beam_skip_error_threshold_"))
    {
        private_nh_.param("beam_skip_error_threshold_", beam_skip_error_threshold_);
    }
    else
    {
        private_nh_.param("beam_skip_error_threshold", beam_skip_error_threshold_, 0.9);
    }

    private_nh_.param("laser_z_hit", z_hit_, 0.95);
    private_nh_.param("laser_z_short", z_short_, 0.1);
    private_nh_.param("laser_z_max", z_max_, 0.05);
    private_nh_.param("laser_z_rand", z_rand_, 0.05);
    private_nh_.param("laser_sigma_hit", sigma_hit_, 0.2);
    private_nh_.param("laser_lambda_short", lambda_short_, 0.1);
    private_nh_.param("laser_likelihood_max_dist", laser_likelihood_max_dist_, 2.0);
    std::string tmp_model_type;
    private_nh_.param("laser_model_type", tmp_model_type, std::string("likelihood_field"));		// 如果没有配置，默认使用 【似然域模型】
    if (tmp_model_type == "beam")
        laser_model_type_ = LASER_MODEL_BEAM;
    else if (tmp_model_type == "likelihood_field")
        laser_model_type_ = LASER_MODEL_LIKELIHOOD_FIELD;
    else if (tmp_model_type == "likelihood_field_prob")
    {
        laser_model_type_ = LASER_MODEL_LIKELIHOOD_FIELD_PROB;
    }
    else
    {
        ROS_WARN("Unknown laser model type \"%s\"; defaulting to likelihood_field model",
                 tmp_model_type.c_str());
        laser_model_type_ = LASER_MODEL_LIKELIHOOD_FIELD;
    }

    private_nh_.param("odom_model_type", tmp_model_type, std::string("diff"));
    if (tmp_model_type == "diff")
        odom_model_type_ = ODOM_MODEL_DIFF;
    else if (tmp_model_type == "omni")
        odom_model_type_ = ODOM_MODEL_OMNI;
    else if (tmp_model_type == "diff-corrected")
        odom_model_type_ = ODOM_MODEL_DIFF_CORRECTED;
    else if (tmp_model_type == "omni-corrected")
        odom_model_type_ = ODOM_MODEL_OMNI_CORRECTED;
    else
    {
        ROS_WARN("Unknown odom model type \"%s\"; defaulting to diff model",
                 tmp_model_type.c_str());
        odom_model_type_ = ODOM_MODEL_DIFF;
    }

    private_nh_.param("update_min_d", d_thresh_, 0.2);
    private_nh_.param("update_min_a", a_thresh_, M_PI / 6.0);
    private_nh_.param("odom_frame_id", odom_frame_id_, std::string("odom"));
    private_nh_.param("base_frame_id", base_frame_id_, std::string("base_link"));
    private_nh_.param("global_frame_id", global_frame_id_, std::string("map"));
    private_nh_.param("resample_interval", resample_interval_, 2);
    private_nh_.param("selective_resampling", selective_resampling_, false);
    double tmp_tol;
    private_nh_.param("transform_tolerance", tmp_tol, 0.1);
    private_nh_.param("recovery_alpha_slow", alpha_slow_, 0.001);
    private_nh_.param("recovery_alpha_fast", alpha_fast_, 0.1);
    private_nh_.param("tf_broadcast", tf_broadcast_, true);
    private_nh_.param("force_update_after_initialpose", force_update_after_initialpose_, false);
    private_nh_.param("force_update_after_set_map", force_update_after_set_map_, false);

    // For diagnostics
    private_nh_.param("std_warn_level_x", std_warn_level_x_, 0.2);
    private_nh_.param("std_warn_level_y", std_warn_level_y_, 0.2);
    private_nh_.param("std_warn_level_yaw", std_warn_level_yaw_, 0.1);

    transform_tolerance_.fromSec(tmp_tol);

    {
        double bag_scan_period;
        private_nh_.param("bag_scan_period", bag_scan_period, -1.0);
        bag_scan_period_.fromSec(bag_scan_period);
    }

    odom_frame_id_ = stripSlash(odom_frame_id_);
    base_frame_id_ = stripSlash(base_frame_id_);
    global_frame_id_ = stripSlash(global_frame_id_);

    updatePoseFromServer();

    cloud_pub_interval.fromSec(1.0);
    tfb_.reset(new tf2_ros::TransformBroadcaster());
    tf_.reset(new tf2_ros::Buffer());
    tfl_.reset(new tf2_ros::TransformListener(*tf_));
    // publisher service subsciber 初始化
    pose_pub_ = nh_.advertise<geometry_msgs::PoseWithCovarianceStamped>("amcl_pose", 2, true);
    particlecloud_pub_ = nh_.advertise<geometry_msgs::PoseArray>("particlecloud", 2, true);
    global_loc_srv_ = nh_.advertiseService("global_localization",
                                           &AmclNode::globalLocalizationCallback,
                                           this);
    nomotion_update_srv_ = nh_.advertiseService("request_nomotion_update", &AmclNode::nomotionUpdateCallback, this);
    set_map_srv_ = nh_.advertiseService("set_map", &AmclNode::setMapCallback, this);

    laser_scan_sub_ = new message_filters::Subscriber<sensor_msgs::LaserScan>(nh_, scan_topic_, 100);
    laser_scan_filter_ =
        new tf2_ros::MessageFilter<sensor_msgs::LaserScan>(*laser_scan_sub_,
                                                           *tf_,
                                                           odom_frame_id_,
                                                           100,
                                                           nh_);
    laser_scan_filter_->registerCallback(boost::bind(&AmclNode::laserReceived,
                                                     this, _1));
    initial_pose_sub_ = nh_.subscribe("initialpose", 2, &AmclNode::initialPoseReceived, this);

    if (use_map_topic_)
    {
        map_sub_ = nh_.subscribe("map", 1, &AmclNode::mapReceived, this);
        ROS_INFO("Subscribed to map topic.");
    }
    else
    {
        requestMap();
    }
    m_force_update = false;

    dsrv_ = new dynamic_reconfigure::Server<amcl::AMCLConfig>(ros::NodeHandle("~"));
    dynamic_reconfigure::Server<amcl::AMCLConfig>::CallbackType cb = boost::bind(&AmclNode::reconfigureCB, this, _1, _2);
    dsrv_->setCallback(cb);

    // 15s timer to warn on lack of receipt of laser scans, #5209
    laser_check_interval_ = ros::Duration(15.0);
    check_laser_timer_ = nh_.createTimer(laser_check_interval_,
                                         boost::bind(&AmclNode::checkLaserReceived, this, _1));

    diagnosic_updater_.setHardwareID("None");
    diagnosic_updater_.add("Standard deviation", this, &AmclNode::standardDeviationDiagnostics);
}

#pragma region
void AmclNode::reconfigureCB(AMCLConfig &config, uint32_t level)
{}

void AmclNode::runFromBag(const std::string &in_bag_fn, bool trigger_global_localization)
{}

void AmclNode::savePoseToServer()   // 仅仅在2个地方调用  1.laserReceived  2.中断函数中
{       // AMCLCMT: 把当前的base->map的位姿保存到参数服务器中的   initial_pose_x,y,a    initial_cov_xx,initial_cov_yy,initial_cov_aa中
    // We need to apply the last transform to the latest odom pose to get
    // the latest map pose to store.  We'll take the covariance from
    // last_published_pose.
    tf2::Transform odom_pose_tf2;
    tf2::convert(latest_odom_pose_.pose, odom_pose_tf2);        // 把latest_odom_pose_中的矢量和旋转转换为一个tf2::Transform对象odom_pose_tf2
    tf2::Transform map_pose = latest_tf_.inverse() * odom_pose_tf2;
                // AMCLCMT: 将机器人坐标系变换到map坐标系下，map_pose(tf2::Transform中有这个位姿)中保存了这个位姿
    double yaw = tf2::getYaw(map_pose.getRotation());

    ROS_DEBUG("Saving pose to server. x: %.3f, y: %.3f", map_pose.getOrigin().x(), map_pose.getOrigin().y());

    private_nh_.setParam("initial_pose_x", map_pose.getOrigin().x());		// map_pose 是base@map的坐标，将其传入参数服务器的initial_pose_x,initial_pose_y initial_pose_z中；
    private_nh_.setParam("initial_pose_y", map_pose.getOrigin().y());
    private_nh_.setParam("initial_pose_a", yaw);
    private_nh_.setParam("initial_cov_xx",
                         last_published_pose.pose.covariance[6 * 0 + 0]);		
    private_nh_.setParam("initial_cov_yy",
                         last_published_pose.pose.covariance[6 * 1 + 1]);
    private_nh_.setParam("initial_cov_aa",
                         last_published_pose.pose.covariance[6 * 5 + 5]);
}

void AmclNode::updatePoseFromServer()
                // 这里是初始化位姿的前身，还处于设置均值和协方差的阶段：   
                // 主要是根据配置文件加载到参数服务器，然后将参数服务器中的配置项，设置到程序变量init_pose_[3],init_cov_[3] 之中
{
    init_pose_[0] = 0.0;
    init_pose_[1] = 0.0;
    init_pose_[2] = 0.0;
    init_cov_[0] = 0.5 * 0.5;
    init_cov_[1] = 0.5 * 0.5;
    init_cov_[2] = (M_PI / 12.0) * (M_PI / 12.0);
    // Check for NAN on input from param server, #5239
    double tmp_pos;
    private_nh_.param("initial_pose_x", tmp_pos, init_pose_[0]);
    if (!std::isnan(tmp_pos))
        init_pose_[0] = tmp_pos;
    else
        ROS_WARN("ignoring NAN in initial pose X position");
    private_nh_.param("initial_pose_y", tmp_pos, init_pose_[1]);
    if (!std::isnan(tmp_pos))
        init_pose_[1] = tmp_pos;
    else
        ROS_WARN("ignoring NAN in initial pose Y position");
    private_nh_.param("initial_pose_a", tmp_pos, init_pose_[2]);
    if (!std::isnan(tmp_pos))
        init_pose_[2] = tmp_pos;
    else
        ROS_WARN("ignoring NAN in initial pose Yaw");

    private_nh_.param("initial_cov_xx", tmp_pos, init_cov_[0]);
    if (!std::isnan(tmp_pos))
        init_cov_[0] = tmp_pos;
    else
        ROS_WARN("ignoring NAN in initial covariance XX");
    private_nh_.param("initial_cov_yy", tmp_pos, init_cov_[1]);
    if (!std::isnan(tmp_pos))
        init_cov_[1] = tmp_pos;
    else
        ROS_WARN("ignoring NAN in initial covariance YY");
    private_nh_.param("initial_cov_aa", tmp_pos, init_cov_[2]);
    if (!std::isnan(tmp_pos))
        init_cov_[2] = tmp_pos;
    else
        ROS_WARN("ignoring NAN in initial covariance AA");
}

void AmclNode::checkLaserReceived(const ros::TimerEvent &event)
{
}

void AmclNode::requestMap()
{
    boost::recursive_mutex::scoped_lock ml(configuration_mutex_);

    // get map via RPC
    nav_msgs::GetMap::Request req;
    nav_msgs::GetMap::Response resp;
    ROS_INFO("Requesting the map...");
    while (!ros::service::call("static_map", req, resp))
    {
        ROS_WARN("Request for map failed; trying again...");
        ros::Duration d(0.5);
        d.sleep();
    }
    handleMapMessage(resp.map);
}

void AmclNode::mapReceived(const nav_msgs::OccupancyGridConstPtr &msg)
{
    if (first_map_only_ && first_map_received_)
    {
        return;
    }

    handleMapMessage(*msg);

    first_map_received_ = true;
}
#pragma endregion

/* ***************** ***************************************************

    它以地图数据消息的引用为输入参数
    构建了地图对象，里程计对象和雷达传感器对象，
    并完成了滤波器的初始化工作。 
    在函数的一开始对配置信号量加锁，并输出一些关于地图尺寸和分辨率的日志。
 ****************************************************/
// WHOLELINE:   地图处理 hadnldeMapMessage()
void AmclNode::handleMapMessage(const nav_msgs::OccupancyGrid &msg)

{
    boost::recursive_mutex::scoped_lock cfl(configuration_mutex_);

    ROS_INFO("Received a %d X %d map @ %.3f m/pix\n",
             msg.info.width,
             msg.info.height,
             msg.info.resolution);

    if (msg.header.frame_id != global_frame_id_)
        ROS_WARN("Frame_id of map received:'%s' doesn't match global_frame_id:'%s'. This could cause issues with reading published topics",
                 msg.header.frame_id.c_str(),
                 global_frame_id_.c_str());

    freeMapDependentMemory();
    // Clear queued laser objects because they hold pointers to the existing
    // map, #5202.
    lasers_.clear();
    lasers_update_.clear();
    frame_to_laser_.clear();

    map_ = convertMap(msg);     		// AMCLCMT 01. 创建地图对象amcl自己的地图格式

#if NEW_UNIFORM_SAMPLING
    // Index of free space
    free_space_indices.resize(0);           // free_space_indices 先把容器清空
    for (int i = 0; i < map_->size_x; i++)
        for (int j = 0; j < map_->size_y; j++)
            if (map_->cells[MAP_INDEX(map_, i, j)].occ_state == -1)     // 如果栅格状态=-1(空闲)，则将这个栅格的位置(i,j)保存到free_space_indices
                free_space_indices.push_back(std::make_pair(i, j));     // AMCLCMT 02. free_space_indices 存储空闲栅格的位置
#endif
    // Create the particle filter		// AMCLCMT 03. pf_alloc() 返回创建的pf_
    pf_ = pf_alloc(min_particles_, max_particles_,
                   alpha_slow_, alpha_fast_,
                   (pf_init_model_fn_t)AmclNode::uniformPoseGenerator,
                   (void *)map_);
    pf_set_selective_resampling(pf_, selective_resampling_);
    pf_->pop_err = pf_err_;
    pf_->pop_z = pf_z_;

    // Initialize the filter  
    updatePoseFromServer();				// AMCLCMT 04. 将参数服务器中的初始位姿，更新到变量init_pose_[3]中，调用一次，重新更新有一次
    pf_vector_t pf_init_pose_mean = pf_vector_zero();
    pf_init_pose_mean.v[0] = init_pose_[0];
    pf_init_pose_mean.v[1] = init_pose_[1];
    pf_init_pose_mean.v[2] = init_pose_[2];
    pf_matrix_t pf_init_pose_cov = pf_matrix_zero();
    pf_init_pose_cov.m[0][0] = init_cov_[0];
    pf_init_pose_cov.m[1][1] = init_cov_[1];
    pf_init_pose_cov.m[2][2] = init_cov_[2];
    pf_init(pf_, pf_init_pose_mean, pf_init_pose_cov);		// AMCLCMT 05. 根据当前机器人的初始位姿，初始化pf
    pf_init_ = false;

    // Instantiate the sensor objects
    // Odometry  // AMCLCMT: 构建里程计对象并配置里程计模型。
    delete odom_;
    odom_ = new AMCLOdom();
    ROS_ASSERT(odom_);
    odom_->SetModel(odom_model_type_, alpha1_, alpha2_, alpha3_, alpha4_, alpha5_);

    // Laser    // AMCLCMT: 构建雷达传感器对象，并根据运行参数laser_model_type_构建不同模型的雷达。
    delete laser_;
    laser_ = new AMCLLaser(max_beams_, map_);
    ROS_ASSERT(laser_);
    if (laser_model_type_ == LASER_MODEL_BEAM)
        laser_->SetModelBeam(z_hit_, z_short_, z_max_, z_rand_,
                             sigma_hit_, lambda_short_, 0.0);
    else if (laser_model_type_ == LASER_MODEL_LIKELIHOOD_FIELD_PROB)
    {
        ROS_INFO("Initializing likelihood field model; this can take some time on large maps...");
        laser_->SetModelLikelihoodFieldProb(z_hit_, z_rand_, sigma_hit_,
                                            laser_likelihood_max_dist_,
                                            do_beamskip_, beam_skip_distance_,
                                            beam_skip_threshold_, beam_skip_error_threshold_);
        ROS_INFO("Done initializing likelihood field model.");
    }
    else
    {
        ROS_INFO("Initializing likelihood field model; this can take some time on large maps...");
        laser_->SetModelLikelihoodField(z_hit_, z_rand_, sigma_hit_,
                                        laser_likelihood_max_dist_);
        ROS_INFO("Done initializing likelihood field model.");
    }

    // In case the initial pose message arrived before the first map,
    // try to apply the initial pose now that the map has arrived.
    applyInitialPose();
}

void AmclNode::freeMapDependentMemory()
{
    if (map_ != NULL)
    {
        map_free(map_);
        map_ = NULL;
    }
    if (pf_ != NULL)
    {
        pf_free(pf_);
        pf_ = NULL;
    }
    delete odom_;
    odom_ = NULL;
    delete laser_;
    laser_ = NULL;
}

/**
 * Convert an OccupancyGrid map message into the internal
 * representation. This allocates a map_t and returns it.
 */
map_t *
AmclNode::convertMap(const nav_msgs::OccupancyGrid &map_msg)
{
    map_t *map = map_alloc();
    ROS_ASSERT(map);

    map->size_x = map_msg.info.width;		// mapserver发布的width 是int类型，size_x也是int类型
    map->size_y = map_msg.info.height;
    map->scale = map_msg.info.resolution;
    map->origin_x = map_msg.info.origin.position.x + (map->size_x / 2) * map->scale;    // 标准map_msg是地图的左下角在map中的坐标为origin，此处转换为地图的中心在世界坐标系中的坐标
    map->origin_y = map_msg.info.origin.position.y + (map->size_y / 2) * map->scale;
    // Convert to player format
    map->cells = (map_cell_t *)malloc(sizeof(map_cell_t) * map->size_x * map->size_y);
    ROS_ASSERT(map->cells);
    for (int i = 0; i < map->size_x * map->size_y; i++)
    {
        if (map_msg.data[i] == 0)
            map->cells[i].occ_state = -1;
        else if (map_msg.data[i] == 100)
            map->cells[i].occ_state = +1;
        else
            map->cells[i].occ_state = 0;
    }

    return map;
}

AmclNode::~AmclNode()
{
    delete dsrv_;
    freeMapDependentMemory();
    delete laser_scan_filter_;
    delete laser_scan_sub_;
    // TODO: delete everything allocated in constructor
}

bool AmclNode::getOdomPose(geometry_msgs::PoseStamped &odom_pose,
                           double &x, double &y, double &yaw,
                           const ros::Time &t, const std::string &f) 
{
    // Get the robot's pose
    geometry_msgs::PoseStamped ident;
    ident.header.frame_id = stripSlash(f);
    ident.header.stamp = t;
    tf2::toMsg(tf2::Transform::getIdentity(), ident.pose); // ident是有个poseStamped消息，用单位变换，给pose的旋转和平移赋值
    try
    {
        this->tf_->transform(ident, odom_pose, odom_frame_id_);//将f坐标系中的点ident，变换到odom坐标系下，将这个新的坐标保存到odom_pose中
        // odom_pose.header.frame_id ="odom",  此处是将base_link(f坐标系)中的原点（没有旋转）转换到odom坐标系下
    }
    catch (const tf2::TransformException &e)
    {
        ROS_WARN("Failed to compute odom pose, skipping scan (%s)", e.what());
        return false;
    }
    x = odom_pose.pose.position.x;
    y = odom_pose.pose.position.y;
    yaw = tf2::getYaw(odom_pose.pose.orientation);

    return true;
}

pf_vector_t
AmclNode::uniformPoseGenerator(void *arg)       // 随机返回地图中的某个空闲的栅格
{
    map_t *map = (map_t *)arg;                  // 把这个指针，转换为map_t指针；
#if NEW_UNIFORM_SAMPLING
    unsigned int rand_index = drand48() * free_space_indices.size();        //drand48()随机返回一个[0,1]之间的浮点数，  意思是：随机找一个索引；
    std::pair<int, int> free_point = free_space_indices[rand_index];        // 取出该随机索引 处 的栅格位置
    pf_vector_t p;
    p.v[0] = MAP_WXGX(map, free_point.first);
    p.v[1] = MAP_WYGY(map, free_point.second);
    p.v[2] = drand48() * 2 * M_PI - M_PI;
#else
    double min_x, max_x, min_y, max_y;

    min_x = (map->size_x * map->scale) / 2.0 - map->origin_x;
    max_x = (map->size_x * map->scale) / 2.0 + map->origin_x;
    min_y = (map->size_y * map->scale) / 2.0 - map->origin_y;
    max_y = (map->size_y * map->scale) / 2.0 + map->origin_y;

    pf_vector_t p;

    ROS_DEBUG("Generating new uniform sample");
    for (;;)
    {
        p.v[0] = min_x + drand48() * (max_x - min_x);
        p.v[1] = min_y + drand48() * (max_y - min_y);
        p.v[2] = drand48() * 2 * M_PI - M_PI;
        // Check that it's a free cell
        int i, j;
        i = MAP_GXWX(map, p.v[0]);
        j = MAP_GYWY(map, p.v[1]);
        if (MAP_VALID(map, i, j) && (map->cells[MAP_INDEX(map, i, j)].occ_state == -1))
            break;
    }
#endif
    return p;
}

bool AmclNode::globalLocalizationCallback(std_srvs::Empty::Request &req,
                                          std_srvs::Empty::Response &res)
{
    if (map_ == NULL)
    {
        return true;
    }
    boost::recursive_mutex::scoped_lock gl(configuration_mutex_);
    ROS_INFO("Initializing with uniform distribution");
    pf_init_model(pf_, (pf_init_model_fn_t)AmclNode::uniformPoseGenerator,
                  (void *)map_);
    ROS_INFO("Global initialisation done!");
    pf_init_ = false;
    return true;
}

// force nomotion updates (amcl updating without requiring motion)
bool AmclNode::nomotionUpdateCallback(std_srvs::Empty::Request &req,
                                      std_srvs::Empty::Response &res)
{
    m_force_update = true;
    // ROS_INFO("Requesting no-motion update");
    return true;
}

bool AmclNode::setMapCallback(nav_msgs::SetMap::Request &req,
                              nav_msgs::SetMap::Response &res)
{
    handleMapMessage(req.map);
    handleInitialPoseMessage(req.initial_pose);
    if (force_update_after_set_map_)
    {
        m_force_update = true;
    }
    res.success = true;
    return true;
}
// WHOLELINE:  laserReceived 
void AmclNode::laserReceived(const sensor_msgs::LaserScanConstPtr &laser_scan) // laserscan订阅的回调
{
    std::string laser_scan_frame_id = stripSlash(laser_scan->header.frame_id);		// laser_scan_frame_id = laser_link
    last_laser_received_ts_ = ros::Time::now(); // 当前时间，赋值给  last_laser_received_ts
    if (map_ == NULL)
    {
        return;
    }
    boost::recursive_mutex::scoped_lock lr(configuration_mutex_);
    int laser_index = -1;
/**************************************************************************************************************
 * step1 : 判断当前激光帧数据的雷达是否已经存在与 frame_to_laser_中，获取激光雷达在map中的序号，保存在laser_index中
 **************************************************************************************************************/
    // Do we have the base->base_laser Tx yet?
    if (frame_to_laser_.find(laser_scan_frame_id) == frame_to_laser_.end())		  // frame_to_laser_ map中没有laser_link的坐标系 的分支
    {
        ROS_DEBUG("Setting up laser %d (frame_id=%s)\n", (int)frame_to_laser_.size(), laser_scan_frame_id.c_str());
        lasers_.push_back(new AMCLLaser(*laser_));			// AMCLCMT: 问题：此处frame_to_laser容器为空，将一个*laser对象，放入到lasers_容器中，为何，用途是什么？
        lasers_update_.push_back(true);						// AMCLCMT: lasers_<AMCLLaser>  lasers_update_<bool>
        laser_index = frame_to_laser_.size();				// 如果只有一个雷达的话，则laser_index = 0 

        geometry_msgs::PoseStamped ident;						
        ident.header.frame_id = laser_scan_frame_id;
        ident.header.stamp = ros::Time();
        tf2::toMsg(tf2::Transform::getIdentity(), ident.pose);		// ident 表示的是laser_link坐标系下的坐标原点；

        geometry_msgs::PoseStamped laser_pose;				// AMCLCMT： laser_pose 保存laser的原点在base坐标系下的位姿
        try
        {
            this->tf_->transform(ident, laser_pose, base_frame_id_);		
        }	
        catch (const tf2::TransformException &e)
        {
            ROS_ERROR("Couldn't transform from %s to %s, "
                      "even though the message notifier is in use",
                      laser_scan_frame_id.c_str(),
                      base_frame_id_.c_str());
            return;
        }

        pf_vector_t laser_pose_v;
        laser_pose_v.v[0] = 89741.pose.position.x;
        laser_pose_v.v[1] = laser_pose.pose.position.y;
        // laser mounting angle gets computed later -> set to 0 here!
        laser_pose_v.v[2] = 0;
        lasers_[laser_index]->SetLaserPose(laser_pose_v);
        ROS_DEBUG("Received laser's pose wrt robot: %.3f %.3f %.3f",
                  laser_pose_v.v[0],
                  laser_pose_v.v[1],
                  laser_pose_v.v[2]);

        frame_to_laser_[laser_scan_frame_id] = laser_index;		// 第一个map插入 	["laser_link", 0]
    }
    else
    {
        // we have the laser pose, retrieve laser index
        laser_index = frame_to_laser_[laser_scan_frame_id];
    }

/**************************************************************************************************************
 * step2 : 获取当前base相对于odom的位姿 
 **************************************************************************************************************/
    // Where was the robot when this scan was taken?
    pf_vector_t pose;       // AMCLCMT pose，保存的是 base在odom中的x,y,yaw
    if (!getOdomPose(latest_odom_pose_, pose.v[0], pose.v[1], pose.v[2],       
                     laser_scan->header.stamp, base_frame_id_))  
                     // getOdomPose() 将base在odom中的位姿，保存在 latest_odom_pose_ 中，所以"latest_odom_pose_"保存的就是base->odom的坐标变换；
		// AMCLCMT: latest_odom_pose_ 保存的是base_link坐标系的原点在odom坐标系中的坐标，也就是机器人在odom中的坐标
		// AMCLCMT: pose 向量保存了base@odom中的x,y,yaw 三个平面位姿

    {
        ROS_ERROR("Couldn't determine robot's pose associated with laser scan");
        return;
    }

    pf_vector_t delta = pf_vector_zero();

/**************************************************************************************************************
 * step3 : 围绕pf_init_是否已经初始化完成，操作不同分支
 * 	3.1 pf_init = true的时候，判断里程计变化量是否需要更新传感器数据 lasers_update_[i] 是否设置为 true
 * 	3.2 pf_init = true 且lasers_update_[laser_index] = 1，里程计更新
 **************************************************************************************************************/
	// 该分支在(重新预处理之后)第二次收到激光雷达数据的时候，就会执行
	// 记录base在odom坐标系中的位姿偏差，保存在delta中
    if (pf_init_)
    {
        // Compute change in pose
        // delta = pf_vector_coord_sub(pose, pf_odom_pose_);
        delta.v[0] = pose.v[0] - pf_odom_pose_.v[0];
        delta.v[1] = pose.v[1] - pf_odom_pose_.v[1];
        delta.v[2] = angle_diff(pose.v[2], pf_odom_pose_.v[2]);

        // See if we should update the filter
        bool update = fabs(delta.v[0]) > d_thresh_ ||
                      fabs(delta.v[1]) > d_thresh_ ||
                      fabs(delta.v[2]) > a_thresh_;		// 如果位姿偏差的delta（三个分量中，有1个超过了变化阈值，就标记需要更新update = true）
        update = update || m_force_update;
        m_force_update = false;				// 此处说明：m_force_update 在其他位置的设置，仅仅一次生效，比如通过服务设置为"无移动时强制更新"=true
											// 此处又重置为false
        // Set the laser update flags		
        if (update)	// AMCLCMT 如果delta已经达到更新的阈值了，则update设置为true， 如果update=true，则所有雷达数据should update 	lasers_[i]
            for (unsigned int i = 0; i < lasers_update_.size(); i++)
                lasers_update_[i] = true;
    }

    bool force_publication = false;
	// pf_init_ = false时，   
	// 三种情况下，会设置pf_init_ = false （1） handlemapmessage()处理了地图消息之后；（2）调用全局定位服务之后；（3）初始位姿设置之后；
	// 所以该分支只在上述三种情况后才执行；
    if (!pf_init_)	
    {
        // Pose at last filter update
        pf_odom_pose_ = pose;

        // Filter is now initialized  如果不进行上述的三种预处理，此后pf_init_就一直为true了
        pf_init_ = true;		

        // Should update sensor data
        for (unsigned int i = 0; i < lasers_update_.size(); i++)
            lasers_update_[i] = true;		// 所有雷达数据的更新全部设置为true，// AMCLCMT 表示全部雷达数据需要更新  should update 

        force_publication = true;

        resample_count_ = 0;
    }
    // If the robot has moved, update the filter
    else if (pf_init_ && lasers_update_[laser_index])		// 
    {
        // printf("pose\n");
        // pf_vector_fprintf(pose, stdout, "%.3f");

        AMCLOdomData odata;
        odata.pose = pose;		// pose是当前帧激光 所确定的base在odom中的位姿；
        // HACK
        // Modify the delta in the action data so the filter gets
        // updated correctly
        odata.delta = delta;	// delta是pose-pf_odom_pose ：是当前帧减去上一帧是的位姿，也就是上一帧时刻后的位移增量

        // Use the action data to update the filter
        odom_->UpdateAction(pf_, (AMCLSensorData *)&odata);

        // Pose at last filter update
        // this->pf_odom_pose = pose;
    }
/**************************************************************************************************************
 * step4 : 传感器更新
 **************************************************************************************************************/
    bool resampled = false;
    // If the robot has moved, update the filter			// AMCLCMT:  运动后更新滤波器
    if (lasers_update_[laser_index])  		// 该激光雷达对应的容器中，更新了
    {
        AMCLLaserData ldata;
        ldata.sensor = lasers_[laser_index];
        ldata.range_count = laser_scan->ranges.size();

        // To account for lasers that are mounted upside-down, we determine the
        // min, max, and increment angles of the laser in the base frame.
        //
        // Construct min and max angles of laser, in the base_link frame.
        tf2::Quaternion q;
        q.setRPY(0.0, 0.0, laser_scan->angle_min);		// q是最小角度的四元数(-3.1415926 对应的q)
        geometry_msgs::QuaternionStamped min_q, inc_q;
        min_q.header.stamp = laser_scan->header.stamp;
        min_q.header.frame_id = stripSlash(laser_scan->header.frame_id);		// min_q 的坐标系是laser
        tf2::convert(q, min_q.quaternion);				// min_q 是最小角度的四元数

        q.setRPY(0.0, 0.0, laser_scan->angle_min + laser_scan->angle_increment);	
        inc_q.header = min_q.header;					// inc_q 的坐标系是laser
        tf2::convert(q, inc_q.quaternion);				// inc_q 是最小角度的四元数+ 1分辨率
        try
        {
            tf_->transform(min_q, min_q, base_frame_id_);		// 转换到base坐标系（laser坐标系的原点在base坐标系中的坐标）
            tf_->transform(inc_q, inc_q, base_frame_id_);		
        }
        catch (const tf2::TransformException &e)
        {
            ROS_WARN("Unable to transform min/max laser angles into base frame: %s",
                     e.what());
            return;
        }

        double angle_min = tf2::getYaw(min_q.quaternion);
        double angle_increment = tf2::getYaw(inc_q.quaternion) - angle_min;

        // wrapping angle to [-pi .. pi]
        angle_increment = fmod(angle_increment + 5 * M_PI, 2 * M_PI) - M_PI;

        ROS_DEBUG("Laser %d angles in base frame: min: %.3f inc: %.3f", laser_index, angle_min, angle_increment);

        // Apply range min/max thresholds, if the user supplied them
        if (laser_max_range_ > 0.0)
            ldata.range_max = std::min(laser_scan->range_max, (float)laser_max_range_);		// ldata.range_max 设置激光距离的最大值：为laser_scan中最值或者是配置值，中偏小的那个
        else
            ldata.range_max = laser_scan->range_max;
        double range_min;
        if (laser_min_range_ > 0.0)
            range_min = std::max(laser_scan->range_min, (float)laser_min_range_);		// range_min :设置的激光测距最小值，配置者和laser_scan中的较小的那个
        else
            range_min = laser_scan->range_min;

        if (ldata.range_max <= 0.0 || range_min < 0.0)
        {
            ROS_ERROR("range_max or range_min from laser is negative! ignore this message.");
            return; // ignore this.
        }

        // The AMCLLaserData destructor will free this memory
        ldata.ranges = new double[ldata.range_count][2];		// AMCLCMT ranges的类型为: double (*ranges)[2];   
        ROS_ASSERT(ldata.ranges);		// ranges[2] 两个元素，第一个表示距离，第二个表示角度；
        for (int i = 0; i < ldata.range_count; i++)
        {
            // amcl doesn't (yet) have a concept of min range.  So we'll map short
            // readings to max range.
            if (laser_scan->ranges[i] <= range_min)
                ldata.ranges[i][0] = ldata.range_max;		// 测量距离小于设置的最小值，则ranges[i][0] = range_max
            else if (laser_scan->ranges[i] > ldata.range_max)		// 测量距离大于设置的最大值，则ranges[i][0] = double的max()
                ldata.ranges[i][0] = std::numeric_limits<decltype(ldata.range_max)>::max();	// decltype() 根据变量进行类型推断； std::numeric_limits<TYPE>::max(); 获取某种类型的最大值
            else
                ldata.ranges[i][0] = laser_scan->ranges[i];		// 满足距离要求
            // Compute bearing
            ldata.ranges[i][1] = angle_min +
                                 (i * angle_increment);		// ranges[i][1] 角度
        }

        lasers_[laser_index]->UpdateSensor(pf_, (AMCLSensorData *)&ldata);

        lasers_update_[laser_index] = false;			// AMCLCMT 对激光扫描数据执行完UpdateSensor之后，就表名该雷达数据已经更新了，所以此处lasers_update[laser_index] 设置为false
														// 就表示不再更新了
        pf_odom_pose_ = pose;

        // Resample the particles
        if (!(++resample_count_ % resample_interval_))
        {
            pf_update_resample(pf_);
            resampled = true;
        }

        pf_sample_set_t *set = pf_->sets + pf_->current_set;
        ROS_DEBUG("Num samples: %d\n", set->sample_count);

        // Publish the resulting cloud
        // TODO: set maximum rate for publishing
        if (!m_force_update)
        {
            geometry_msgs::PoseArray cloud_msg;
            cloud_msg.header.stamp = ros::Time::now();
            cloud_msg.header.frame_id = global_frame_id_;
            cloud_msg.poses.resize(set->sample_count);
            for (int i = 0; i < set->sample_count; i++)
            {
                cloud_msg.poses[i].position.x = set->samples[i].pose.v[0];
                cloud_msg.poses[i].position.y = set->samples[i].pose.v[1];
                cloud_msg.poses[i].position.z = 0;
                tf2::Quaternion q;
                q.setRPY(0, 0, set->samples[i].pose.v[2]);
                tf2::convert(q, cloud_msg.poses[i].orientation);
            }
            particlecloud_pub_.publish(cloud_msg);
        }
    }

    if (resampled || force_publication)
    {
        // Read out the current hypotheses
        double max_weight = 0.0;
        int max_weight_hyp = -1;
        std::vector<amcl_hyp_t> hyps;       // amcl_hyp_t 容器
        hyps.resize(pf_->sets[pf_->current_set].cluster_count);
        for (int hyp_count = 0;
             hyp_count < pf_->sets[pf_->current_set].cluster_count; hyp_count++)
        {
            double weight;
            pf_vector_t pose_mean;
            pf_matrix_t pose_cov;
            if (!pf_get_cluster_stats(pf_, hyp_count, &weight, &pose_mean, &pose_cov))
            {
                ROS_ERROR("Couldn't get stats on cluster %d", hyp_count);
                break;
            }

            hyps[hyp_count].weight = weight;
            hyps[hyp_count].pf_pose_mean = pose_mean;
            hyps[hyp_count].pf_pose_cov = pose_cov;

            if (hyps[hyp_count].weight > max_weight)
            {
                max_weight = hyps[hyp_count].weight;
                max_weight_hyp = hyp_count;     // 最大权重对应的粒子索引
            }
        }

        if (max_weight > 0.0)       // 如果最大权重能够找出来
        {
            ROS_DEBUG("Max weight pose: %.3f %.3f %.3f",
                      hyps[max_weight_hyp].pf_pose_mean.v[0],
                      hyps[max_weight_hyp].pf_pose_mean.v[1],
                      hyps[max_weight_hyp].pf_pose_mean.v[2]);

            /*
               puts("");
               pf_matrix_fprintf(hyps[max_weight_hyp].pf_pose_cov, stdout, "%6.3f");
               puts("");
             */

            geometry_msgs::PoseWithCovarianceStamped p;     // 带协方差的位姿p，用于保存最大权重粒子的位姿均值？
            // Fill in the header
            p.header.frame_id = global_frame_id_;
            p.header.stamp = laser_scan->header.stamp;
            // Copy in the pose
            p.pose.pose.position.x = hyps[max_weight_hyp].pf_pose_mean.v[0];
            p.pose.pose.position.y = hyps[max_weight_hyp].pf_pose_mean.v[1];

            tf2::Quaternion q;
            q.setRPY(0, 0, hyps[max_weight_hyp].pf_pose_mean.v[2]);
            tf2::convert(q, p.pose.pose.orientation);
            // Copy in the covariance, converting from 3-D to 6-D
            pf_sample_set_t *set = pf_->sets + pf_->current_set;
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    // Report the overall filter covariance, rather than the
                    // covariance for the highest-weight cluster
                    // p.covariance[6*i+j] = hyps[max_weight_hyp].pf_pose_cov.m[i][j];
                    p.pose.covariance[6 * i + j] = set->cov.m[i][j];
                }
            }
            // Report the overall filter covariance, rather than the
            // covariance for the highest-weight cluster
            // p.covariance[6*5+5] = hyps[max_weight_hyp].pf_pose_cov.m[2][2];
            p.pose.covariance[6 * 5 + 5] = set->cov.m[2][2];

            /*
               printf("cov:\n");
               for(int i=0; i<6; i++)
               {
               for(int j=0; j<6; j++)
               printf("%6.3f ", p.covariance[6*i+j]);
               puts("");
               }
             */

            pose_pub_.publish(p);			// AMCLCMT : 发布最大权重的集群的pose统计值
            last_published_pose = p;

            ROS_DEBUG("New pose: %6.3f %6.3f %6.3f",
                      hyps[max_weight_hyp].pf_pose_mean.v[0],
                      hyps[max_weight_hyp].pf_pose_mean.v[1],
                      hyps[max_weight_hyp].pf_pose_mean.v[2]);

            // subtracting base to odom from map to base and send map to odom instead
            geometry_msgs::PoseStamped odom_to_map;
            try
            {
                tf2::Quaternion q;
                q.setRPY(0, 0, hyps[max_weight_hyp].pf_pose_mean.v[2]);
				// AMCLCMT：tmp_tf是base_link在global map下的坐标，即base-->map
                tf2::Transform tmp_tf(q, tf2::Vector3(hyps[max_weight_hyp].pf_pose_mean.v[0],
                                                      hyps[max_weight_hyp].pf_pose_mean.v[1],
                                                      0.0));        // 

                geometry_msgs::PoseStamped tmp_tf_stamped;
                tmp_tf_stamped.header.frame_id = base_frame_id_;        // tmp_tf_stamped   坐标系为base
                tmp_tf_stamped.header.stamp = laser_scan->header.stamp;
                tf2::toMsg(tmp_tf.inverse(), tmp_tf_stamped.pose);      // 将tmp_tf变换取inverse(),得到的是map->base_link的变换； 将这个变换的位姿存入tmp_tf_stamped.pose中；

                this->tf_->transform(tmp_tf_stamped, odom_to_map, odom_frame_id_);      // 应该是把map->base的坐标，变换成map->odom的坐标，将这个位姿保存在odom_to_map
            }       // 所以odom_to_map ：保存的是map在odom下的位姿
            catch (const tf2::TransformException &)
            {
                ROS_DEBUG("Failed to subtract base to odom transform");
                return;
            }
		// tf2::Transform latest_tf_;  构造函数中的声明  
            tf2::convert(odom_to_map.pose, latest_tf_);     // latest_tf_ 保存的是map->odom的变换；
            latest_tf_valid_ = true;

            if (tf_broadcast_ == true)			// 如果要发布tf的话
            {
                // We want to send a transform that is good up until a
                // tolerance time so that odom can be used
                ros::Time transform_expiration = (laser_scan->header.stamp +
                                                  transform_tolerance_);
                geometry_msgs::TransformStamped tmp_tf_stamped;
                tmp_tf_stamped.header.frame_id = global_frame_id_;		// map
                tmp_tf_stamped.header.stamp = transform_expiration;
                tmp_tf_stamped.child_frame_id = odom_frame_id_;
				// tmp_tf_stamped 表示 parent:map  child:odom    即 odom->map 的坐标变换
                tf2::convert(latest_tf_.inverse(), tmp_tf_stamped.transform);

                this->tfb_->sendTransform(tmp_tf_stamped);		// 把odom->map 的tf发布出来
                sent_first_transform_ = true;
            }
        }
        else
        {
            ROS_ERROR("No pose!");
        }
    }
    else if (latest_tf_valid_)
    {
        if (tf_broadcast_ == true)
        {
            // Nothing changed, so we'll just republish the last transform, to keep
            // everybody happy.
            ros::Time transform_expiration = (laser_scan->header.stamp +
                                              transform_tolerance_);
            geometry_msgs::TransformStamped tmp_tf_stamped;
            tmp_tf_stamped.header.frame_id = global_frame_id_;
            tmp_tf_stamped.header.stamp = transform_expiration;
            tmp_tf_stamped.child_frame_id = odom_frame_id_;
            tf2::convert(latest_tf_.inverse(), tmp_tf_stamped.transform);		// AMCLCMT ：tmp_tf_stamped这个变换是odom原点在map坐标系的坐标，即odom-->map
            this->tfb_->sendTransform(tmp_tf_stamped);
        }

        // Is it time to save our last pose to the param server
        ros::Time now = ros::Time::now();
        if ((save_pose_period.toSec() > 0.0) &&
            (now - save_pose_last_time) >= save_pose_period)
        {
            this->savePoseToServer();
            save_pose_last_time = now;
        }
    }

    diagnosic_updater_.update();
}

void AmclNode::initialPoseReceived(const geometry_msgs::PoseWithCovarianceStampedConstPtr &msg)
{
    handleInitialPoseMessage(*msg);
    if (force_update_after_initialpose_)
    {
        m_force_update = true;
    }
}

void AmclNode::handleInitialPoseMessage(const geometry_msgs::PoseWithCovarianceStamped &msg)
{
    boost::recursive_mutex::scoped_lock prl(configuration_mutex_);
    if (msg.header.frame_id == "")
    {
        // This should be removed at some point
        ROS_WARN("Received initial pose with empty frame_id.  You should always supply a frame_id.");
    }
    // We only accept initial pose estimates in the global frame, #5148.
    else if (stripSlash(msg.header.frame_id) != global_frame_id_)
    {
        ROS_WARN("Ignoring initial pose in frame \"%s\"; initial poses must be in the global frame, \"%s\"",
                 stripSlash(msg.header.frame_id).c_str(),
                 global_frame_id_.c_str());
        return;
    }

    // In case the client sent us a pose estimate in the past, integrate the
    // intervening odometric change.
    geometry_msgs::TransformStamped tx_odom;
    try
    {
        // wait a little for the latest tf to become available		// 调用了 lookuptransform的高级接口，含有6个参数
        tx_odom = tf_->lookupTransform(base_frame_id_, msg.header.stamp,		// 现在时间，与收到的激光雷达的时间戳是有区别的，此处是获取，在当前时间与时间戳上相差的这一点时间内，base_link的变化
                                       base_frame_id_, ros::Time::now(),		// 此行是子坐标系，所以表示的是，当前时刻在时间戳时的坐标系下的坐标变换关系
                                       odom_frame_id_, ros::Duration(0.5));		// 把odom坐标系作为固定坐标系
    }	// tx_odom 现在的base相对于时间戳里的base 的坐标关系
    catch (const tf2::TransformException &e)
    {
        // If we've never sent a transform, then this is normal, because the
        // global_frame_id_ frame doesn't exist.  We only care about in-time
        // transformation for on-the-move pose-setting, so ignoring this
        // startup condition doesn't really cost us anything.
        if (sent_first_transform_)
            ROS_WARN("Failed to transform initial pose in time (%s)", e.what());
        tf2::convert(tf2::Transform::getIdentity(), tx_odom.transform);
    }

    tf2::Transform tx_odom_tf2;
    tf2::convert(tx_odom.transform, tx_odom_tf2);		// tx_odom_tf2: 是base_link之间（两个时间戳之间的base_link是有位姿上的差异的）的相对关系
    tf2::Transform pose_old, pose_new;
    tf2::convert(msg.pose.pose, pose_old);				// msg是相对于map坐标系的一个初始位姿
    pose_new = pose_old * tx_odom_tf2;

    // Transform into the global frame

    ROS_INFO("Setting pose (%.6f): %.3f %.3f %.3f",
             ros::Time::now().toSec(),
             pose_new.getOrigin().x(),
             pose_new.getOrigin().y(),
             tf2::getYaw(pose_new.getRotation()));
    // Re-initialize the filter
    pf_vector_t pf_init_pose_mean = pf_vector_zero();
    pf_init_pose_mean.v[0] = pose_new.getOrigin().x();
    pf_init_pose_mean.v[1] = pose_new.getOrigin().y();
    pf_init_pose_mean.v[2] = tf2::getYaw(pose_new.getRotation());
    pf_matrix_t pf_init_pose_cov = pf_matrix_zero();
    // Copy in the covariance, converting from 6-D to 3-D
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            pf_init_pose_cov.m[i][j] = msg.pose.covariance[6 * i + j];
        }
    }
    pf_init_pose_cov.m[2][2] = msg.pose.covariance[6 * 5 + 5];

    delete initial_pose_hyp_;
    initial_pose_hyp_ = new amcl_hyp_t();
    initial_pose_hyp_->pf_pose_mean = pf_init_pose_mean;
    initial_pose_hyp_->pf_pose_cov = pf_init_pose_cov;
    applyInitialPose();
}

/**
 * If initial_pose_hyp_ and map_ are both non-null, apply the initial
 * pose to the particle filter state.  initial_pose_hyp_ is deleted
 * and set to NULL after it is used.
 */
void AmclNode::applyInitialPose()		// 使用获得的初始位姿，初始化粒子滤波器pf_
{
    boost::recursive_mutex::scoped_lock cfl(configuration_mutex_);
    if (initial_pose_hyp_ != NULL && map_ != NULL)      // initial_pose_hyp_非空和map_地图对象也有的情况下，才执行后续代码
    {
        pf_init(pf_, initial_pose_hyp_->pf_pose_mean, initial_pose_hyp_->pf_pose_cov);      // 用initial_pose_hyp_ 变量初始化粒子滤波器
        pf_init_ = false;       // ??? 表示还没有初始化？？？

        delete initial_pose_hyp_;
        initial_pose_hyp_ = NULL;
    }
}

void AmclNode::standardDeviationDiagnostics(diagnostic_updater::DiagnosticStatusWrapper &diagnostic_status)
{
    double std_x = sqrt(last_published_pose.pose.covariance[6 * 0 + 0]);        // pose.pose.covariance[36]是协方差矩阵，写成了1维形式，std_x 是x的标准差
    double std_y = sqrt(last_published_pose.pose.covariance[6 * 1 + 1]);        // std_y 是y的标准差
    double std_yaw = sqrt(last_published_pose.pose.covariance[6 * 5 + 5]);       // std_yaw 偏航角yaw的标准差

    diagnostic_status.add("std_x", std_x);
    diagnostic_status.add("std_y", std_y);
    diagnostic_status.add("std_yaw", std_yaw);
    diagnostic_status.add("std_warn_level_x", std_warn_level_x_);
    diagnostic_status.add("std_warn_level_y", std_warn_level_y_);
    diagnostic_status.add("std_warn_level_yaw", std_warn_level_yaw_);

    if (std_x > std_warn_level_x_ || std_y > std_warn_level_y_ || std_yaw > std_warn_level_yaw_)
    {
        diagnostic_status.summary(diagnostic_msgs::DiagnosticStatus::WARN, "Too large");
    }
    else
    {
        diagnostic_status.summary(diagnostic_msgs::DiagnosticStatus::OK, "OK");
    }
}
