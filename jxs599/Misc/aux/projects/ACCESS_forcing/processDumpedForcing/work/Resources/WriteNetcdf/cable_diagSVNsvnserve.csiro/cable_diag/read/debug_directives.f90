!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!### uncomment the model in which you are running
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
#define UM 
!#define offline
!
#define comp_data_directive
#define plot_data_directive
!
#ifdef UM
#define DEFr_1 8
#define DEFi_d 8
#endif
!
#ifdef offline 
#define DEFr_1 kind(1.1)
#define DEFi_d kind(1)
#endif
!
!
#define gs kind(1.1)
#define gi kind(1)
!
!#define testxs
!
!#define manual_extrema
#ifdef manual_extrema
#define fymin -1.0
#define fymax 1.5
#endif


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
