# Rotation Matrix Reconstruction in Quaternions with an Optimal Observer: A Navigation Algorithm Based on Nonlinear Complementary Observers

## Ariel Iporre Rivas∗and Mauricio Amestegui Moreno ́ 

```
Electronic Engineering - Faculty of Engineering - Universidad Mayor de San Andres ́
La Paz, Bolivia
```
Abstract—This work proposes a variation of an earlier nav-
igation algorithm in the literature composed by SO(3) comple-
mentary observers. Such modification means the inclusion of a
quaternion optimal observer to determine the rotation matrix
instead using the original vectorial reconstruction approach. In
particular, this proposal involves an experimental comparison
between the original and modified method. The results show
about 40% increase for estimation quality while 21% more
complexity using this new approach.
Index Terms—Inertial navigation, Kalman filters, Quaternions,
INS, Navigation Systems, Optimal Observers, Nonlinear Filters

REFERENCES
[1] A. Lukyanov, S. Dodds, and J. Vittek, “Observer-based attitude control
in the sliding mode,”WIT Transactions on The Built Environment,
vol. 22, 1970.
[2] S. Nicosia and P. Tomei, “Nonlinear observer and output feedback
attitude control of spacecraft,”IEEE Transactions on Aerospace and
Electronic Systems, vol. 28, no. 4, pp. 970–977, 1992.
[3] J. F. Vasconcelos, C. Silvestre, and P. Oliveira, “A nonlinear gps/imu
based observer for rigid body attitude and position estimation,” in
Decision and Control, 2008. CDC 2008. 47th IEEE Conference on.
IEEE, 2008, pp. 1255–1260.
[4] G. G. Scandaroli and P. Morin, “Nonlinear filter design for pose and
imu bias estimation,” inRobotics and Automation (ICRA), 2011 IEEE
International Conference on. IEEE, 2011, pp. 4524–4530.
[5] R. Mahony, T. Hamel, and J.-M. Pflimlin, “Nonlinear complementary
filters on the special orthogonal group,”IEEE Transactions on automatic
control, vol. 53, no. 5, pp. 1203–1218, 2008.
[6] R. E. Kalman, “A new approach to linear filtering and prediction
problems,”Transactions of the ASME–Journal of Basic Engineering,
vol. 82, no. Series D, pp. 35–45, 1960.
[7] S. Schmidt, G. Smith, and L. Mcgee, “Application of statistical filter
theory to the optimal estimation of position and velocity on board a
circumlunar vehicle,” Ames Research Center, Tech. Rep., 1962.
[8] V. Bistrovs and A. Kluga, “Adaptive extended kalman filter for aided
inertial navigation system,”Elektronika ir Elektrotechnika, vol. 122,
no. 6, pp. 37–40, 2012.
[9] F. A. Faruqi and K. J. Turner, “Extended kalman filter synthesis
for integrated global positioning/inertial navigation systems,”Applied
mathematics and computation, vol. 115, no. 2, pp. 213–227, 2000.
[10] F. Woyano, S. Lee, and S. Park, “Evaluation and comparison of
performance analysis of indoor inertial navigation system based on foot
mounted IMU,” in2016 18th International Conference on Advanced
Communication Technology (ICACT), Jan 2016, pp. 792–798.
[11] J. L. Marins, X. Yun, E. R. Bachmann, R. B. McGhee, and M. J. Zyda,
“An extended kalman filter for quaternion-based orientation estimation
using marg sensors,” inIntelligent Robots and Systems, 2001. Proceed-
ings. 2001 IEEE/RSJ International Conference on, vol. 4. IEEE, 2001,
pp. 2003–2011.
[12] D. Luenberger, “Observers for multivariable systems,”IEEE Transac-
tions on Automatic Control, vol. 11, no. 2, pp. 190–197, 1966.
[13] F. E. Thau, “Observing the state of non-linear dynamic systems,”
International Journal of Control, vol. 17, no. 3, pp. 471–479, 1973.
[Online]. Available: http://dx.doi.org/10.1080/
[14] J. Thienel and R. M. Sanner, “A coupled nonlinear spacecraft attitude
controller and observer with an unknown constant gyro bias and gyro
noise,”IEEE Transactions on Automatic Control, vol. 48, no. 11, pp.
2011–2015, 2003.
[15] M.-D. Hua, “Attitude observers for accelerated rigid bodies based on
gps and ins measurements,” inDecision and Control, 2009 held jointly
with the 2009 28th Chinese Control Conference. CDC/CCC 2009.
Proceedings of the 48th IEEE Conference on. IEEE, 2009, pp. 8071–
8076.
[16] S. Bonnabel, P. Martin, and P. Rouchon, “Symmetry-preserving ob-
servers,”IEEE Transactions on Automatic Control, vol. 53, no. 11, pp.
2514–2526, 2008.
[17] J. Sola, “Quaternion kinematics for the error-state kf,” 2015.
[18] F. L. Lewis, D. Vrabie, and V. L. Syrmos,Optimal control. John Wiley
& Sons, 2012.
