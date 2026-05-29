#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX 105
#define EPS 1e-12

void printMenu(){
    printf("      CHUONG TRINH TINH TOAN       \n");
    printf("====================================\n");
    printf("1. Giai he phuong trinh theo phuong phap Krame\n");
    printf("2. Giai he phuong trinh theo phuong phap Gauss\n");
    printf("3. Giai he phuong trinh theo phuong phap lap Gauss-Seidel\n");
    printf("4. Giai he phuong trinh theo phuong phap Giam du\n");
    printf("5. Tim gia tri rieng va vector rieng\n");
    printf("6. Noi suy Lagrange\n");
    printf("7. Noi suy Lagrange voi cac moc cach deu\n");
    printf("8. Noi suy Ayken\n");
    printf("9. Noi suy Ayken loai 2\n");
    printf("10. Noi suy Newton\n");
    printf("11. Noi suy Hecmit\n");
    printf("12.Phuong phap binh phuong be nhat\n");
    printf("0. Thoat\n");
    printf("====================================\n");
    printf("Chon: ");
}

//Krame
void Krame(int n, double A[MAX][MAX], double B[MAX][MAX], double X[MAX]){
    int i, j;
    printf("\n========================================\n");
    printf("             MA TRAN HE SO              \n");
    printf("========================================\n");
    for (i=1; i<=n; i++){
        printf("|");
        for (j=1; j<=n+1; j++){
            printf("%10.2f", A[i][j]);
        }
        printf("  |\n");
    }
    int cottThuK, iTam, jTam, kTam, pTam, jdoidong;
    double tmp;
    for (cottThuK=0; cottThuK<=n; cottThuK++){
        for (i=1; i<=n; i++){
            for (j=1; j<=n; j++){
                if (j==cottThuK){
                    B[i][j] = A[i][n+1];
                } else{
                    B[i][j] = A[i][j];
                }
            }
        }
        int doidaudinhthuc = 0;
        double dinhthuc = 1.0;
        int isZero = 0;
        for (iTam=1; iTam<n; iTam++){
            if (B[iTam][iTam]==0){
                for (kTam = iTam+1; kTam<=n; kTam++){
                    if (B[kTam][iTam]!=0){
                        for (jdoidong=1; jdoidong<=n; jdoidong++){
                            tmp = B[iTam][jdoidong];
                            B[iTam][jdoidong] = B[kTam][jdoidong];
                            B[kTam][jdoidong] = tmp;
                        }
                        doidaudinhthuc++;
                        break;
                    }
                }
            }
            if (B[iTam][iTam]==0){
                isZero = 1;
                break;
            }
            for (jTam=iTam+1; jTam<=n; jTam++){
                double m = -B[jTam][iTam] / B[iTam][iTam];
                for (pTam=1; pTam<=n; pTam++){
                    B[jTam][pTam] = B[iTam][pTam] * m + B[jTam][pTam];
                }
            }
        }
        if (isZero==1){
            X[cottThuK] = 0;
        } else{
            for (int tTam=1; tTam<=n; tTam++){
                dinhthuc *= B[tTam][tTam];
            }
            X[cottThuK] = dinhthuc * pow(-1, doidaudinhthuc);
        }

        printf("Det(A%d) = %.2f\n", cottThuK, X[cottThuK]);

        if (cottThuK==0 && X[0]==0){
            break;
        }
    }
    printf("\n");
    if (X[0]==0){
        printf(" He phuong trinh vo nghiem hoac vo so nghiem\n");
    } else{
        printf("      NGHIEM CUA HE PHUONG TRINH        \n");
        for (int kp=1; kp<=n; kp++){
            printf("x%d = %.4f\n", kp, X[kp] / X[0]);
        }
    }
    printf("\n");
}

// Gauss
void doidong(int n, double A[MAX][MAX], int p, int q){
    for (int k=1; k<=n+1; k++){
        double temp = A[p][k];
        A[p][k] = A[q][k];
        A[q][k] = temp;
    }
}

void Gauss(int n, double A[MAX][MAX], double x[MAX]){
    int i, j, k;
    int kt_nghiem = 1; //1 là có nghiệm duy nhất, 0 là vô nghiệm/vô số nghiệm
    printf("             MA TRAN HE SO              \n");
    printf("========================================\n");
    for (i=1; i<=n; i++){
        printf("|");
        for (j=1; j<=n+1; j++){
            printf("%10.2f", A[i][j]);
        }
        printf("  |\n");
    }
    for (i=1; i<=n-1; i++){
        if (A[i][i]==0){
            int kt = 0;
            for (j=i+1; j<=n; j++){
                if (A[j][i] != 0){
                    doidong(n, A, i, j);
                    kt = 1;
                    break;
                }
            }
            if (!kt){
                kt_nghiem = 0;
                break;
            }
        }

        if (kt_nghiem==1){
            for (j=i+1; j<=n; j++){
                double m = -A[j][i] / A[i][i];
                for (k=i; k<=n+1; k++){
                    A[j][k] = A[j][k] + A[i][k] * m;
                }
            }
        }
    }
    if (A[n][n]==0) kt_nghiem = 0;
    if (kt_nghiem==1){
        for (i=n; i>=1; i--){
            x[i] = A[i][n+1];
            for (j=i+1; j<=n; j++){
                x[i] = x[i] - A[i][j] * x[j];
            }
            x[i] = x[i] / A[i][i];
        }
        printf("      NGHIEM CUA HE PHUONG TRINH        \n");
        printf("========================================\n");
        for (i=1; i<=n; i++){
            printf("x%d = %.4f\n", i, x[i]);
        }
    } else{
        printf(" He phuong trinh vo nghiem hoac vo so nghiem\n");
    }
    printf("\n");
}

// Gauss-Seidel
void Gauss_Seidel(int n, double A[MAX][MAX]){
    double x[MAX], y[MAX];
    int i, j, t;
    double s;

    int iter = 0;
    int max_iter = 10000;

    for (i=1; i <= n; i++){
        x[i] = 0.0;
    }
    printf("             MA TRAN HE SO              \n");
    printf("========================================\n");
    for (i=1; i<=n; i++){
        printf("|");
        for (j=1; j<=n+1; j++){
            printf("%10.2f", A[i][j]);
        }
        printf("  |\n");
    }
    for (i=1; i<=n; i++){
        if (A[i][i]==0){
            printf("\nLOI! A[%d][%d] = 0. Khong the chia cho 0!\n\n", i, i);
            return;
        }
    }
    do{
        t = 0;
        iter++;
        for (i=1; i<=n; i++){
            s = 0;
            for (j=1; j<=n; j++){
                if (j!=i){
                    s = s + A[i][j] * x[j];
                }
            }

            y[i] = (A[i][n + 1] - s) / A[i][i];

            if (isnan(y[i]) || isinf(y[i])){
                printf("\nMa tran phan ky, tran so tai vong lap %d!\n\n", iter);
                return;
            }

            if (fabs(y[i]-x[i]) >= EPS){
                t = 1;
            }
            x[i] = y[i];
        }
        if (iter > max_iter){
            printf("\nKhong hoi tu sau %d vong lap!\n\n", max_iter);
            return;
        }
    } while (t);
    printf("\nDa hoi tu thanh cong sau %d phep lap.\n", iter);
    printf("      NGHIEM CUA HE PHUONG TRINH        \n");
    printf("========================================\n");
    for (i=1; i<=n; i++){
        printf("x%d = %.6f\n", i, x[i]);
    }
    printf("\n");
}

// Giảm dư
void Giam_Du(int n, double A[MAX][MAX], double x[MAX]){
    int i, j, k, t, iter = 0;
    int max_iter = 10000;
    double a_tmp[MAX][MAX];
    double r[MAX];
    double d, imax, temp;

    for (i=1; i<=n; i++){
        x[i] = 0.0;
    }
    printf("             MA TRAN HE SO              \n");
    printf("========================================\n");
    for (i=1; i<=n; i++){
        printf("|");
        for (j=1; j<=n+1; j++){
            printf("%10.2f", A[i][j]);
        }
        printf("  |\n");
    }
    for (i=1; i<=n; i++){
        if (A[i][i]==0){
            printf("\nLOI! Phan tu A[%d][%d] = 0. Khong the chia cho 0!\n\n", i, i);
            return;
        }
        for (j=1; j<=n+1; j++){
            a_tmp[i][j] = A[i][j];
        }
    }
    for (i=1; i<=n; i++){
        temp = a_tmp[i][i];
        for (j=1; j <=n+1; j++){
            a_tmp[i][j] = a_tmp[i][j] / temp;
        }
    }
    for (i=1; i<=n; i++){
        r[i] = a_tmp[i][n+1];
        for (j=1; j<=n; j++){
            r[i] = r[i] - a_tmp[i][j] * x[j];
        }
    }
    do{
        t = 0;
        iter++;

        imax = fabs(r[1]);
        k = 1;
        for (i=2; i<=n; i++){
            if (imax < fabs(r[i])){
                imax = fabs(r[i]);
                k = i;
            }
        }
        x[k] = x[k] + r[k];
        d = r[k];
        for (i=1; i<=n; i++){
            r[i] = r[i] - a_tmp[i][k] * d;
            if (fabs(r[i]) >= EPS){
                t = 1;
            }
        }
        if (isnan(x[k]) || isinf(x[k])){
            printf("\nMa tran phan ky, tran so tai vong lap %d!\n\n", iter);
            return;
        }
        if (iter > max_iter){
            printf("\nKhong hoi tu sau %d vong lap!\n\n", max_iter);
            return;
        }
    } while (t);
    printf("\nDa hoi tu thanh cong sau %d phep lap.\n", iter);
    printf("      NGHIEM CUA HE PHUONG TRINH        \n");
    printf("========================================\n");
    for (i=1; i<=n; i++){
        printf("x%d = %.6f\n", i, x[i]);
    }
    printf("\n");
}

//GIA TRI RIENG VA VECTO RIENG
void nhan2mt(int n, double A[MAX][MAX], double B[MAX][MAX], double C[MAX][MAX]){
    for (int i=1; i<=n; i++){
        for (int j=1; j<=n; j++){
            C[i][j]=0;
            for (int k=1; k<=n; k++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void danhilepski(int n, double A[MAX][MAX], double C[MAX][MAX]){
    double m[MAX][MAX], m1[MAX][MAX], B[MAX][MAX];
    for (int i=1; i<=n; i++){
        for (int j=1; j<=n; j++){
            C[i][j] = (i==j) ? 1.0 : 0.0; // Ma tran don vi
        }
    }
    for (int k=n-1; k>=1; k--){
        for (int i=1; i<=n; i++){
            for (int j=1; j<=n; j++){
                if (i!=k){
                    m[i][j] = (i==j) ? 1.0 : 0.0;
                    m1[i][j] = (i==j) ? 1.0 : 0.0;
                } else {
                    m1[i][j] = A[k + 1][j];
                    if (j==k) m[i][j] = 1.0 / A[k + 1][k];
                    else m[i][j] = -A[k + 1][j] / A[k + 1][k];
                }
            }
        }
        nhan2mt(n, A, m, B);
        nhan2mt(n, m1, B, A);
        nhan2mt(n, C, m, B);
        for (int i=1; i<=n; i++)
            for (int j=1; j<=n; j++)
                C[i][j] = B[i][j];
    }
}

double tinh_f_lambda(double A[MAX][MAX], int n, double x_ct){
    double p = pow(x_ct, n);
    for (int i=1; i<=n; i++){
        p -= A[1][i] * pow(x_ct, n-i);
    }
    return p;
}

void tinh_x(double C[MAX][MAX], int n, double lambda, double x_vec[MAX][MAX], int *idx){
    double y[MAX], x[MAX];
    for (int i=1; i<=n; i++) y[i] = pow(lambda, n-i);
    for (int i=1; i<=n; i++){
        x[i]=0;
        for (int j=1; j<=n; j++){
            x[i] += C[i][j] * y[j];
        }
    }
    for (int i=1; i<=n; i++) x_vec[*idx][i] = x[i];
    (*idx)++;
}

void Danilevsky_Full(int n, double A[MAX][MAX]){
    double a_tmp[MAX][MAX], C[MAX][MAX];
    double x_vec[MAX][MAX];
    int idx = 1;
    printf("             MA TRAN BAN DAU            \n");
    printf("========================================\n");
    for (int i=1; i<=n; i++){
        printf("|");
        for (int j=1; j<=n; j++){
            a_tmp[i][j] = A[i][j];
            printf("%10.2f", a_tmp[i][j]);
        }
        printf("  |\n");
    }

    danhilepski(n, a_tmp, C);
    printf("             MA TRAN P \n");
    printf("========================================\n");
    for(int i=1; i<=n; i++){
        printf("|");
        for(int j=1; j<=n; j++) printf("%10.2f", a_tmp[i][j]);
        printf("  |\n");
    }
    double max_p = 0;
    for (int i=1; i<=n; i++){
        if (fabs(a_tmp[1][i]) > max_p) max_p = fabs(a_tmp[1][i]);
    }
    double R = 1.0 + max_p;
    double b_step = 0.1;
    printf("            CAC GIA TRI RIENG           \n");
    printf("========================================\n");

    int dem = 0;
    for (double x = -R; x <= R; x += b_step){
        double x1 = x;
        double x2 = x + b_step;
        double f1 = tinh_f_lambda(a_tmp, n, x1);
        double f2 = tinh_f_lambda(a_tmp, n, x2);

        if (f1 * f2 < 0){
            double o;
            double a1 = x1, b2 = x2;
            while ((b2 - a1) >= EPS){
                o = (a1 + b2) / 2.0;
                if (tinh_f_lambda(a_tmp, n, o) == 0.0) break;
                else if (tinh_f_lambda(a_tmp, n, a1) * tinh_f_lambda(a_tmp, n, o) < 0) b2 = o;
                else a1 = o;
            }
            dem++;
            double lambda = (a1 + b2) / 2.0;
            printf("  Lambda_%d = %-8.4f\n", dem, lambda);
            tinh_x(C, n, lambda, x_vec, &idx);
        }
    }
    if (dem == 0){
        printf("Khong tim thay gia tri rieng thuc nao trong khoang [-%.2f, %.2f]\n", R, R);
        return;
    }
    printf("            CAC VECTOR RIENG            \n");
    printf("========================================\n");
    for (int i=1; i<idx; i++){
        printf("  v%d = (", i);
        for (int j=1; j<=n; j++){
            printf("%7.3f", x_vec[i][j]);
            if (j<n) printf(", ");
        }
        printf(" )^T\n");
    }
    printf("\n");
}
//LAGRANGE
void inDaThucLagrange(int n, double x[MAX], double y[MAX]){
    printf("       BIEU THUC DA THUC LAGRANGE       \n");
    printf("========================================\n");
    printf("P(x) = \n");

    int isFirst = 1;
    for (int i=0; i<n; i++){
        if (y[i] == 0) continue;
        if (!isFirst){
            printf("  +  ");
        }
        isFirst = 0;
        double mau_so = 1.0;
        for (int j=0; j<n; j++){
            if (i!=j) mau_so *= (x[i] - x[j]);
        }
        printf("(%.2f / %.2f)", y[i], mau_so);
        for (int j=0; j<n; j++){
            if (i!=j){
                if (x[j] > 0) printf(" * (x - %.2f)", x[j]);
                else if (x[j] < 0) printf(" * (x + %.2f)", fabs(x[j]));
                else printf(" * x");
            }
        }
        printf("\n");
    }
}
void tinhToanLagrange(int n, double x[MAX], double y[MAX], double x_tinh, double *kqCach1, double *kqCach2){
    int i, j;
    *kqCach1 = 0.0;
    *kqCach2 = 0.0;
// TÍNH CÁCH 1
    for (i=0; i<n; i++){
        double pi = 1.0;
        for (j=0; j<n; j++){
            if (j!=i){
                pi = pi * (x_tinh - x[j]) / (x[i] - x[j]);
            }
        }
        *kqCach1 += y[i] * pi;
    }
// TÍNH CÁCH 2
    double W = 1.0;
    for (i=0; i<n; i++){
        W *= (x_tinh - x[i]);
    }

    for (i=0; i<n; i++){
        double Wi = 1.0;
        for (j=0; j<n; j++){
            if (j!=i){
                Wi *= (x[i] - x[j]);
            }
        }
        *kqCach2 += y[i] / ((x_tinh - x[i]) * Wi);
    }
    *kqCach2 *= W;
}

void noiSuyLagrange(int n, double x[MAX], double y[MAX]){
    int i;
    double x_tinh, kq1, kq2;
    printf("          NOI SUY LAGRANGE              \n");
    printf("========================================\n");
    printf("            CAC MOC NOI SUY             \n");
    printf("========================================\n");
    for (i=0; i<n; i++){
        printf("x = %-10.2f y = %-10.2f\n", x[i], y[i]);
    }
    inDaThucLagrange(n, x, y);
    printf("\nNhap gia tri x can tinh f(x): ");
    scanf("%lf", &x_tinh);
    int isNode = -1;
    double minX = x[0], maxX = x[0];
    for (i=0; i<n; i++){
        if (x[i] < minX) minX = x[i];
        if (x[i] > maxX) maxX = x[i];
        if (fabs(x_tinh - x[i]) < EPS) isNode = i;
    }
    if (isNode!=-1){
        printf("\nx can tinh trung voi moc noi suy x[%d].\n", isNode);
        printf("f(%.2f) = %.6f\n\n", x_tinh, y[isNode]);
    }
    else {
        if (x_tinh < minX || x_tinh > maxX) {
            printf("\nx can tinh nam ngoai khoang mien noi suy.\n");
        }
        tinhToanLagrange(n, x, y, x_tinh, &kq1, &kq2);
        printf("\n================ KET QUA ===============\n");
        printf(" Cach 1: f(%.2f) = %.6f\n", x_tinh, kq1);
        printf(" Cach 2: f(%.2f) = %.6f\n", x_tinh, kq2);
    }
}

//Lagrange voi cac moc cach deu
void nhan(double d[MAX], int bac, double k, double kq[MAX]){
    for (int i=0; i<=bac+1; i++){
        kq[i] = 0;
    }
    for (int i=0; i<=bac; i++){
        kq[i] += -k * d[i];
        kq[i+1] += d[i];
    }
}

int kiemTraCachDeu(double x[MAX], int n, double *h){
    if (n<=1){
        *h = 0;
        return 1;
    }
    *h = x[1] - x[0];
    if (fabs(*h) < EPS) return 0;

    for (int i=2; i<n; i++){
        double l = x[i] - x[i-1];
        if (fabs(l - *h) > 1e-9) return 0;
    }
    return 1;
}

int lagrangeCachDeu(double x[MAX], double y[MAX], int n, double p[MAX]){
    double h;
    if (!kiemTraCachDeu(x, n, &h)){
        return 0;
    }
    for (int i=0; i<n; i++) {
        p[i] = 0;
    }
    for (int i=0; i<n; i++) {
        double g[MAX] = {1}; // Mảng g lưu đa thức tử số
        int bac = 0;
        double mau = 1;
        for (int j=0; j<n; j++){
            if (i!=j){
                mau *= (x[i] - x[j]);
                double tam[MAX];

                // Đã sửa 'd' thành 'g'
                nhan(g, bac, x[j], tam);
                bac++;
                for (int k=0; k<=bac; k++){
                    // Đã sửa 'l' thành 'g'
                    g[k] = tam[k];
                }
            }
        }
        double hs = y[i] / mau;
        for (int k=0; k<=bac; k++) {
            p[k] += g[k] * hs;
        }
    }
    return 1;
}

void inVaTinhDaThuc(int n, double x[MAX], double y[MAX], double p[MAX]){
    if (!lagrangeCachDeu(x, y, n, p)){
        printf("\nCac moc x khong cach deu.\n");
        return;
    }
    printf("       DA THUC NOI SUY LAGRANGE         \n");
    printf("========================================\n");
    printf("P(x) = ");

    int isFirst = 1;
    for (int i=n-1; i>=0; i--){
        if (fabs(p[i]) < EPS) continue;
        if (p[i] > 0 && !isFirst) printf(" + ");
        else if (p[i] < 0 && !isFirst) printf(" - ");
        else if (p[i] < 0 && isFirst) printf("-");

        printf("%.4f", fabs(p[i]));

        if (i>0){
            printf("x");
            if (i>1) printf("^%d", i);
        }
        isFirst = 0;
    }
    printf("\n");

    double x_tinh;
    printf("Nhap gia tri x can tinh P(x): ");
    scanf("%lf", &x_tinh);

    double ket_qua = 0.0;
    for (int i=n-1; i>=0; i--){
        ket_qua = ket_qua * x_tinh + p[i];
    }

    printf("\n=> KET QUA: P(%.4f) = %.6f\n\n", x_tinh, ket_qua);
}

//AYKEN
void noiSuyAyken(int n, double x[MAX], double y[MAX]){
    double c;
    double diff[MAX][MAX];
    double d[MAX];
    int i, j;
    printf("\nNhap gia tri x can tinh f(x): ");
    scanf("%lf", &c);
    printf("                              BANG NOI SUY AYKEN                          \n");
    printf("==========================================================================\n");
    printf("%-5s | %-8s | ", "i", "y_i");
    for (j=0; j<n; j++){
        printf("x_i - x_%-2d | ", j);
    }
    printf("%-10s\n", "D_i");
    printf("--------------------------------------------------------------------------\n");

    for (i=0; i<n; i++){
        d[i] = 1.0;
        printf("x_%-3d | %-8.3f | ", i, y[i]);

        for (j = 0; j < n; j++) {
            if (i == j) {
                diff[i][j] = c - x[i];
            } else {
                diff[i][j] = x[i] - x[j];
            }
            printf("%-10.3f | ", diff[i][j]);
            d[i] *= diff[i][j];
        }
        printf("%-10.3f\n", d[i]);
    }
    printf("==========================================================================\n");

    double w = 1.0;
    for (i=0; i<n; i++){
        w *= (c - x[i]);
    }
    printf("\n=> W(%.3f) = %.6f\n", c, w);

    double s = 0.0;
    for (i=0; i<n; i++){
        if (fabs(c - x[i])<EPS){
            printf("\nDiem can tinh trung voi moc x[%d].\n", i);
            printf("\n================ KET QUA ===============\n");
            printf(" f(%.3f) = %.6f\n", c, y[i]);
            return;
        }
        if (fabs(d[i]) < EPS){
            printf("\nCac moc x bi trung lap (D_i = 0)\n");
            return;
        }
        s += y[i] / d[i];
    }
    printf("Tong S  = %.6f\n", s);
    double kq = w * s;
    printf("\n================ KET QUA ===============\n");
    printf(" f(%.3f) = W * S = %.6f\n", c, kq);
}

//AYKEN LOẠI 2
void noiSuyAykenDang2(int n, double x[MAX], double y[MAX]){
    double x_tinh;
    double P[MAX][MAX];
    int i, j;
    printf("\nNhap gia tri x can tinh f(x): ");
    scanf("%lf", &x_tinh);

    for (i=0; i<n; i++){
        P[i][0] = y[i];
    }
    for (j=1; j<n; j++){
        for (i=j; i<n; i++){
            double tu_so = P[j-1][j-1] * (x[i] - x_tinh) - P[i][j-1] * (x[j-1] - x_tinh);
            double mau_so = x[i] - x[j-1];
            if (fabs(mau_so)<EPS){
                printf("\nCo 2 moc x[%d] va x[%d] trung nhau.\n", i, j-1);
                return;
            }
            P[i][j] = tu_so / mau_so;
        }
    }
    printf("                              BANG NOI SUY AYKEN (DANG 2)                                 \n");
    printf("==========================================================================================\n");
    printf("%-5s | %-8s | ", "x_i", "y_i");
    for (j=1; j<n; j++){
        if (j == 1) printf("%-8s | ", "L_0i");
        else if (j == 2) printf("%-8s | ", "L_01i");
        else if (j == 3) printf("%-8s | ", "L_012i");
        else if (j == 4) printf("%-8s | ", "L_0123i");
        else printf("L_0..%-2di | ", j);
    }
    printf("%-8s\n", "x_i - x");
    printf("------------------------------------------------------------------------------------------\n");

    for (i=0; i<n; i++){
        printf("%-5.1f | %-8.3f | ", x[i], P[i][0]);

        for (j=1; j<n; j++){
            if (j<=i){
                printf("%-8.4f | ", P[i][j]);
            } else {
                printf("%-8s | ", "");
            }
        }
        printf("%-8.2f\n", x[i] - x_tinh);
    }
    printf("==========================================================================================\n");
    printf("\n================ KET QUA ===============\n");
    printf(" => f(%.3f) = %.6f\n", x_tinh, P[n-1][n-1]);
}

//NỘI SUY NEWTON
void noiSuyNewton(int n, double x[MAX], double y[MAX]){
    double P[MAX][MAX];
    int i, j;
    double x_tinh, h;
    h = x[1] - x[0];
    if (fabs(h) < EPS){
        printf("\nBuoc nhay h = 0. Cac moc x khong duoc trung nhau!\n");
        return;
    }
    for (i=2; i<n; i++){
        if (fabs((x[i]-x[i-1])- h)>1e-9){
            printf("\nCac moc x khong cach deu!\n");
            return;
        }
    }
    for (i=0; i<n; i++){
        P[i][0] = y[i];
    }
    for (j=1; j<n; j++){
        for (i=j; i<n; i++){
            P[i][j] = P[i][j-1] - P[i-1][j-1];
        }
    }
    printf("                              BANG SAI PHAN NEWTON                                 \n");
    printf("==========================================================================================\n");

    printf("%-5s | %-8s | ", "x_i", "f(x_i)");
    for (j = 1; j < n; j++) {
        printf("d^%-2d f   | ", j);
    }
    printf("\n------------------------------------------------------------------------------------------\n");
    for (i=0; i<n; i++){
        printf("%-5.1f | %-8.3f | ", x[i], P[i][0]);
        for (j=1; j<n; j++){
            if (j<=i){
                printf("%-8.3f | ", P[i][j]);
            } else {
                printf("%-8s | ", "");
            }
        }
        printf("\n");
    }
    printf("==========================================================================================\n");
    printf("\nNhap gia tri x can tinh f(x): ");
    scanf("%lf", &x_tinh);

    double kq = P[0][0];
    double tich = 1.0;
    double giai_thua = 1.0;
    double h_pow = 1.0;

    for (int k = 1; k < n; k++) {
        tich *= (x_tinh - x[k-1]);
        giai_thua *= k;
        h_pow *= h;
        kq += (P[k][k] / (giai_thua * h_pow)) * tich;
    }

    printf("\n================ KET QUA ===============\n");
    printf(" => f(%.3f) = %.6f\n", x_tinh, kq);
}

//NỘI SUY HECMIT
void noiSuyHecmit(int n, double x[MAX], double y[MAX], int has_dy[MAX], double dy[MAX]){
    double Z[MAX];
    double Q[MAX][MAX] = {0};
    int m = 0;
    int i, j;
    double x_tinh;

    for (i=0; i<n; i++){
        Z[m] = x[i];
        Q[m][0] = y[i];
        m++;

        if (has_dy[i]==1){
            Z[m] = x[i];
            Q[m][0] = y[i];
            Q[m][1] = dy[i];
            m++;
        }
    }
    for (j=1; j<m; j++){
        for (i=j; i<m; i++){
            if (j==1 && fabs(Z[i] - Z[i-1])<EPS){
                continue;
            }
            Q[i][j] = (Q[i][j-1] - Q[i-1][j-1]) / (Z[i] - Z[i-j]);
        }
    }
    printf("                             BANG SAI PHAN CHIA HECMIT                                  \n");
    printf("==========================================================================================\n");
    printf("%-5s | %-8s | ", "Z_i", "f(Z_i)");
    for (j = 1; j < m; j++) printf("f[..%-2d] | ", j);
    printf("\n------------------------------------------------------------------------------------------\n");

    for (i = 0; i < m; i++) {
        printf("%-5.1f | %-8.3f | ", Z[i], Q[i][0]);
        for (j = 1; j < m; j++) {
            if (j <= i) printf("%-8.3f | ", Q[i][j]);
            else printf("%-8s | ", "");
        }
        printf("\n");
    }
    printf("==========================================================================================\n");
    printf("\n               DA THUC NOI SUY HECMIT                  \n");
    printf("========================================================\n");
    printf("H(x) = %.4f\n", Q[0][0]);
    for (i = 1; i < m; i++) {
        if (fabs(Q[i][i]) < EPS) continue;

        if (Q[i][i] > 0) printf("     + %.4f", Q[i][i]);
        else printf("     - %.4f", fabs(Q[i][i]));

        for (int k = 0; k < i; k++) {
            if (Z[k] == 0) printf(" * x");
            else if (Z[k] > 0) printf(" * (x - %.1f)", Z[k]);
            else printf(" * (x + %.1f)", fabs(Z[k]));
        }
        printf("\n");
    }
    printf("========================================================\n");

    printf("\nNhap gia tri x can tinh H(x): ");
    scanf("%lf", &x_tinh);
    double ket_qua = 0.0;
    for (i = m - 1; i >= 0; i--) {
        double tich = 1.0;
        for (int k = 0; k < i; k++) {
            tich *= (x_tinh - Z[k]);
        }
        ket_qua += Q[i][i] * tich;
    }

    printf("\n================ KET QUA ===============\n");
    printf(" => H(%.3f) = %.6f\n", x_tinh, ket_qua);
    printf("========================================\n\n");
}

//PHƯƠNG PHÁP BÌNH PHƯƠNG BÉ NHẤT
void giaiHePT_Mini(int m, double A[MAX][MAX], double X[MAX]){
    int i, j, k;
    for (k = 1; k <= m - 1; k++) {
        for (i = k + 1; i <= m; i++) {
            double factor = A[i][k] / A[k][k];
            for (j = k; j <= m + 1; j++) {
                A[i][j] -= factor * A[k][j];
            }
        }
    }
    for (i = m; i >= 1; i--) {
        X[i] = A[i][m + 1];
        for (j = i + 1; j <= m; j++) {
            X[i] -= A[i][j] * X[j];
        }
        X[i] /= A[i][i];
    }
}

void binhPhuongBeNhat(int n, double x[MAX], double y[MAX]){
    int choice;
    printf("           PHUONG PHAP BINH PHUONG BE NHAT              \n");
    printf("========================================================\n");
    printf("Chon dang ham thuc nghiem can xap xi:\n");
    printf("  1. Tuyen tinh      : y = a + bx\n");
    printf("  2. Phi tuyen bac 2 : y = a + bx + cx^2\n");
    printf("  3. Ham mu          : y = a * e^(bx)\n");
    printf("  4. Ham luy thua    : y = a * x^b\n");
    printf("========================================================\n");
    printf("Nhap lua chon cua ban (1-4): ");
    scanf("%d", &choice);

    double sum_X = 0, sum_Y = 0, sum_XX = 0, sum_XY = 0;
    double sum_XXX = 0, sum_XXXX = 0, sum_XXY = 0;
    double X_arr[MAX], Y_arr[MAX];
    int i;

    for (i = 0; i < n; i++) {
        if (choice == 1 || choice == 2) {
            X_arr[i] = x[i];
            Y_arr[i] = y[i];
        } else if (choice == 3) {
            if (y[i] <= 0) {
                printf("\nHam mu yeu cau y > 0. Kiem tra lai du lieu!\n\n");
                return;
            }
            X_arr[i] = x[i];
            Y_arr[i] = log(y[i]);
        } else if (choice == 4) {
            if (x[i] <= 0 || y[i] <= 0) {
                printf("\nHam luy thua yeu cau x > 0 va y > 0.\n\n");
                return;
            }
            X_arr[i] = log10(x[i]);
            Y_arr[i] = log10(y[i]);
        }
    }
    for (i = 0; i < n; i++) {
        sum_X += X_arr[i];
        sum_Y += Y_arr[i];
        sum_XX += X_arr[i] * X_arr[i];
        sum_XY += X_arr[i] * Y_arr[i];

        if (choice == 2) {
            sum_XXX += X_arr[i] * X_arr[i] * X_arr[i];
            sum_XXXX += X_arr[i] * X_arr[i] * X_arr[i] * X_arr[i];
            sum_XXY += X_arr[i] * X_arr[i] * Y_arr[i];
        }
    }
    double HPT[MAX][MAX], Nghiem[MAX];

    printf("\n=> Bang cac tong Sigma:\n");
    printf("   Sum_X = %.5f, Sum_Y = %.5f\n", sum_X, sum_Y);
    printf("   Sum_X^2 = %.5f, Sum_XY = %.5f\n", sum_XX, sum_XY);

    if (choice == 2) {
        printf("   Sum_X^3 = %.5f, Sum_X^4 = %.5f, Sum_X^2Y = %.5f\n", sum_XXX, sum_XXXX, sum_XXY);
        HPT[1][1] = n;       HPT[1][2] = sum_X;   HPT[1][3] = sum_XX;   HPT[1][4] = sum_Y;
        HPT[2][1] = sum_X;   HPT[2][2] = sum_XX;  HPT[2][3] = sum_XXX;  HPT[2][4] = sum_XY;
        HPT[3][1] = sum_XX;  HPT[3][2] = sum_XXX; HPT[3][3] = sum_XXXX; HPT[3][4] = sum_XXY;

        giaiHePT_Mini(3, HPT, Nghiem);
        double a = Nghiem[1], b = Nghiem[2], c = Nghiem[3];

        printf("\n================ KET QUA ===============\n");
        printf("He so a = %.5f\n", a);
        printf("He so b = %.5f\n", b);
        printf("He so c = %.5f\n", c);
        printf("\nCong thuc thuc nghiem phu hop nhat:\n");
        printf(" y = %.5f + %.5f*x + %.5f*x^2\n", a, b, c);
        printf("========================================\n\n");

    } else {
        HPT[1][1] = n;     HPT[1][2] = sum_X;  HPT[1][3] = sum_Y;
        HPT[2][1] = sum_X; HPT[2][2] = sum_XX; HPT[2][3] = sum_XY;

        giaiHePT_Mini(2, HPT, Nghiem);
        double A = Nghiem[1], B = Nghiem[2];

        printf("\n================ KET QUA ===============\n");
        if (choice == 1) {
            printf("He so a = %.5f\n", A);
            printf("He so b = %.5f\n", B);
            printf("\nCong thuc thuc nghiem phu hop nhat:\n");
            printf(" y = %.5f + %.5f*x\n", A, B);
        }
        else if (choice == 3) {
            double a = exp(A);
            double b = B;
            printf("He so A (Lna) = %.5f\n", A);
            printf("He so B       = %.5f\n", B);
            printf("=> a = e^A    = %.5f\n", a);
            printf("=> b = B      = %.5f\n", b);
            printf("\nCong thuc thuc nghiem phu hop nhat:\n");
            printf(" => y = %.5f * e^(%.5f * x)\n", a, b);
        }
        else if (choice == 4) {
            double a = pow(10, A);
            double b = B;
            printf("He so A (Lga) = %.5f\n", A);
            printf("He so B       = %.5f\n", B);
            printf("=> a = 10^A   = %.5f\n", a);
            printf("=> b = B      = %.5f\n", b);
            printf("\nCong thuc thuc nghiem phu hop nhat:\n");
            printf(" => y = %.5f * x^(%.5f)\n", a, b);
        }
        printf("========================================\n\n");
    }
}

// ĐỌC FILE HỆ PT
int docFilehept(char filename[MAX], int *n, double A[MAX][MAX]){
    FILE *fp;
    int i, j;

    fp = fopen(filename, "r");
    if (fp==NULL){
        printf("\nKhong the mo file '%s'. Vui long kiem tra lai ten file!\n", filename);
        return 0;
    }
    fscanf(fp, "%d", n);
    for (i=1; i<=*n; i++){
        for (j=1; j<=*n+1; j++){
            fscanf(fp, "%lf", &A[i][j]);
        }
    }
    fclose(fp);
    printf("\nDa doc du lieu tu file thanh cong!\n");
    return 1;
}

// ĐỌC FILE MA TRẬN VUÔNG
int docFileMaTranVuong(char filename[MAX], int *n, double A[MAX][MAX]){
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("\nKhong the mo file '%s'. Vui long kiem tra lai!\n", filename);
        return 0;
    }
    if (fscanf(fp, "%d", n)!=1){
        printf("\nFile rong hoac sai dinh dang o dong dau tien!\n");
        fclose(fp);
        return 0;
    }
    for (int i=1; i<=*n; i++){
        for (int j=1; j<=*n; j++){
            fscanf(fp, "%lf", &A[i][j]);
        }
    }
    fclose(fp);
    printf("\nDa doc du lieu ma tran vuong tu file thanh cong!\n");
    return 1;
}

// ĐỌC FILE LAGRANGE
int docFileLagrange(char filename[MAX], int *n, double x[MAX], double y[MAX]) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("\nKhong the mo file '%s'. Vui long kiem tra lai ten hoac vi tri file!\n", filename);
        return 0;
    }
    if (fscanf(fp, "%d", n) != 1){
        printf("\nFile rong hoac sai dinh dang o dong dau tien!\n");
        fclose(fp);
        return 0;
    }
    if (*n <= 0 || *n >= MAX) {
        printf("\nSo luong moc noi suy (%d) khong hop le hoac qua lon!\n", *n);
        fclose(fp);
        return 0;
    }
    int i;
    for (i=0; i<*n; i++){
        fscanf(fp, "%lf", &x[i]);
    }
    for (i=0; i<*n; i++){
        fscanf(fp, "%lf", &y[i]);
    }
    fclose(fp);
    printf("\nDa doc %d moc noi suy tu file!\n", *n);
    return 1;
}

//ĐỌC FILE HECMIT
int docFileHecmit(char filename[MAX], int *n, double x[MAX], double y[MAX], int has_dy[MAX], double dy[MAX]){
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("\nKhong the mo file '%s'. Vui long kiem tra lai!\n", filename);
        return 0;
    }
    if (fscanf(fp, "%d", n) != 1) return 0;
    for (int i = 0; i < *n; i++) fscanf(fp, "%lf", &x[i]);
    for (int i = 0; i < *n; i++) fscanf(fp, "%lf", &y[i]);
    for (int i = 0; i < *n; i++) fscanf(fp, "%d", &has_dy[i]);
    for (int i = 0; i < *n; i++) fscanf(fp, "%lf", &dy[i]);

    fclose(fp);
    printf("\nDa doc thanh cong du lieu Hecmit tu file!\n");
    return 1;
}



int main(){
    int choose, n;
    double A[MAX][MAX], B[MAX][MAX], X[MAX];
    double x[MAX], y[MAX], p[MAX], dy[MAX];
    char filename[MAX];
    int has_dy[MAX];

    do{
        printMenu();
        scanf("%d",&choose);
        if(choose>=1 && choose<=12){
            printf("Nhap file: ");
            scanf("%s",filename);
        }
         switch(choose){
            case 1:
                if(docFilehept(filename,&n,A)) Krame(n,A,B,X);
                break;
            case 2:
                if(docFilehept(filename,&n,A)) Gauss(n,A,x);
                break;
            case 3:
                if(docFilehept(filename,&n,A)) Gauss_Seidel(n,A);
                break;
            case 4:
                if(docFilehept(filename,&n,A)) Giam_Du(n,A,x);
                break;
            case 5:
                if(docFileMaTranVuong(filename,&n,A)) Danilevsky_Full(n,A);
                break;
            case 6:
                if(docFileLagrange(filename,&n,x,y)) noiSuyLagrange(n,x,y);
                break;
            case 7:
                if(docFileLagrange(filename,&n,x,y)) inVaTinhDaThuc(n, x, y, p);
                break;
            case 8:
                if(docFileLagrange(filename,&n,x,y)) noiSuyAyken(n, x, y);
                break;
            case 9:
                if(docFileLagrange(filename,&n,x,y)) noiSuyAykenDang2(n, x, y);
                break;
            case 10:
                if(docFileLagrange(filename,&n,x,y)) noiSuyNewton(n, x, y);
                break;
            case 11:
                if(docFileHecmit(filename,&n,x,y,has_dy,dy)) noiSuyHecmit(n,x,y,has_dy,dy);
                break;
            case 12:
                if (docFileLagrange(filename, &n, x, y)) binhPhuongBeNhat(n, x, y);
            break;
            case 0:
                printf("\nCam on ban da su dung chuong trinh\n");
                break;
            default:
                printf("\nLua chon khong hop le!\n");
        }

    }while(choose!=0);

    return 0;
}
