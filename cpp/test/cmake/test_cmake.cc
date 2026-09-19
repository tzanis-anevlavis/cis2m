#include <Eigen/Dense>

#include <cis2m/brunovskytransformation.hpp>
#include <cis2m/cis_generator.hpp>
#include <cis2m/hpolyhedron.hpp>

int main() {
    const Eigen::MatrixXd A = Eigen::MatrixXd::Zero(1, 1);
    const Eigen::MatrixXd B = Eigen::MatrixXd::Ones(1, 1);
    const cis2m::BrunovskyTransformation transformation(A, B);

    const Eigen::MatrixXd G = Eigen::MatrixXd::Identity(1, 1);
    const Eigen::VectorXd F = Eigen::VectorXd::Zero(1);
    const cis2m::HPolyhedron polyhedron(G, F);

    return (transformation.IsValid() && polyhedron.IsValid()) ? 0 : 1;
}
