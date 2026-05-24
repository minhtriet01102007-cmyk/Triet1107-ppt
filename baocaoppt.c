#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPS 0.001
#define MAX 100

double tinh_mu(double x[], int mu_x, double y[], int mu_y, int n) {
    double s = 0;
    for (int i = 0; i < n; i++) {
        s += pow(x[i], mu_x) * pow(y[i], mu_y);
    }
    return s;
}

void in_bang_hoi_quy(double x[], double y[], int n, int bac_2) {
    printf("\n+-----+------------+------------+------------+------------+------------+\n");
    printf("|  i  |     x      |     y      |    x^2     |     xy     |   x^2*y    |\n");
    printf("+-----+------------+------------+------------+------------+------------+\n");
    for (int i = 0; i < n; i++) {
        printf("| %3d | %10.4f | %10.4f | %10.4f | %10.4f | %10.4f |\n",
               i + 1, x[i], y[i], x[i]*x[i], x[i]*y[i], x[i]*x[i]*y[i]);
    }
    printf("+-----+------------+------------+------------+------------+------------+\n");
    printf("|Tong | %10.4f | %10.4f | %10.4f | %10.4f | %10.4f |\n",
           tinh_mu(x,1,y,0,n), tinh_mu(x,0,y,1,n), tinh_mu(x,2,y,0,n), tinh_mu(x,1,y,1,n), tinh_mu(x,2,y,1,n));
    printf("+-----+------------+------------+------------+------------+------------+\n");
}

// CRAMER
void Cramer(double A[MAX][MAX], int n) {
    double B[MAX][MAX], X[MAX], dinhthuc;
    int doidaudinhthuc, ketqua = 0;

    for (int cottThuK = 0; cottThuK <= n; cottThuK++) {
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                B[i][j] = (j == cottThuK) ? A[i][n + 1] : A[i][j];
            }
        }

        doidaudinhthuc = 0;
        for (int i = 1; i < n; i++) {
            if (B[i][i] == 0) {
                for (int k = i + 1; k <= n; k++) {
                    if (B[k][i] != 0) {
                        for (int jdoidong = 1; jdoidong <= n; jdoidong++) {
                            double tmp = B[i][jdoidong];
                            B[i][jdoidong] = B[k][jdoidong];
                            B[k][jdoidong] = tmp;
                        }
                        doidaudinhthuc++;
                        break;
                    }
                }
            }
            if (B[i][i] == 0) { ketqua = 10; break; }
            for (int j = i + 1; j <= n; j++) {
                double m = -B[j][i] / B[i][i];
                for (int p = 1; p <= n; p++) B[j][p] += B[i][p] * m;
            }
        }
        if (ketqua == 10) break;

        dinhthuc = 1;
        for (int t = 1; t <= n; t++) dinhthuc *= B[t][t];
        X[cottThuK] = dinhthuc * pow(-1, doidaudinhthuc);
        if (X[0] == 0) break;
    }

    if (ketqua == 10 || X[0] == 0) {
        printf("\nHe phuong trinh vo nghiem hoac vo so nghiem.\n");
    } else {
        printf("\n=== NGHIEM CUA HE(KRAMER) ===\n");
        for (int i = 1; i <= n; i++) printf("x%d = %.4f\n", i, X[i] / X[0]);
    }
}

// GAUSS
void Gauss(double a[MAX][MAX], int n) {
    double x[MAX];
    for (int i = 1; i <= n - 1; i++) {
        for (int j = i + 1; j <= n; j++) {
            if(fabs(a[i][i]) < 1e-9) {
                printf("\nLOI.\n");
                return;
            }
            double m = a[j][i] / a[i][i];
            for (int k = i; k <= n + 1; k++) a[j][k] -= m * a[i][k];
        }
    }
    for (int i = n; i >= 1; i--) {
        double s = a[i][n + 1];
        for (int j = i + 1; j <= n; j++) s -= a[i][j] * x[j];
        x[i] = s / a[i][i];
    }
    printf("\n=== NGHIEM CUA HE(GAUSS) ===\n");
    for (int i = 1; i <= n; i++) printf("x%d = %.6f\n", i, x[i]);
}

// GAUSS-SEIDEL
void GaussSeidel(double a[MAX][MAX], int n) {
    double x[MAX], y[MAX];
    printf("\nNhap nghiem xap xi ban dau cho %d an:\n", n);
    for (int i = 1; i <= n; i++) {
        printf("x%d^(0) = ", i);
        scanf("%lf", &x[i]);
    }

    for(int i=1; i<=n; i++) printf("------------+");
    printf("\n| Lan  |");
    for(int i=1; i<=n; i++) printf("     x%d    |", i);
    printf("\n+------+");
    for(int i=1; i<=n; i++) printf("------------+");

    printf("\n| %4d |", 0);
    for(int i=1; i<=n; i++) printf(" %10.6f |", x[i]);

    int lan_lap = 0, stop;
    do {
        stop = 0;
        lan_lap++;
        for (int i = 1; i <= n; i++) {
            double s = 0;
            for (int j = 1; j <= n; j++) {
                if (i != j) s += a[i][j] * x[j];
            }
            y[i] = (a[i][n + 1] - s) / a[i][i];
            if (fabs(y[i] - x[i]) > EPS) stop = 1;
            x[i] = y[i];
        }

        printf("\n| %4d |", lan_lap);
        for(int i=1; i<=n; i++) printf(" %10.6f |", x[i]);

        if(lan_lap > 500) break;
    } while (stop);

    printf("\n+------+");
    for(int i=1; i<=n; i++) printf("------------+");
    printf("\n=> Nghiem hoi tu sau %d lan lap.\n", lan_lap);
}

// GIẢM DƯ
void GiamDu(double a[MAX][MAX], int n) {
    double x[MAX], r[MAX];
    printf("\nNhap nghiem xap xi ban dau cho %d an:\n", n);
    for (int i = 1; i <= n; i++) {
        printf("x%d^(0) = ", i);
        scanf("%lf", &x[i]);
    }

    for (int i = 1; i <= n; i++) {
        double t = a[i][i];
        for (int j = 1; j <= n + 1; j++) a[i][j] /= t;
    }

    for (int i = 1; i <= n; i++) {
        r[i] = a[i][n + 1];
        for (int j = 1; j <= n; j++) r[i] -= a[i][j] * x[j];
    }

    for(int i=1; i<=n; i++) printf("------------+");
    for(int i=1; i<=n; i++) printf("------------+");
    printf("\n| Vong |");
    for(int i=1; i<=n; i++) printf("     x%d    |", i);
    for(int i=1; i<=n; i++) printf("     r%d    |", i);
    printf("\n+------+");
    for(int i=1; i<=n; i++) printf("------------+");
    for(int i=1; i<=n; i++) printf("------------+");

    int lan_lap = 0, stop;
    do {
        printf("\n| %4d |", lan_lap);
        for(int i=1; i<=n; i++) printf(" %10.6f |", x[i]);
        for(int i=1; i<=n; i++) printf(" %10.6f |", r[i]);

        stop = 0;
        int k = 1;
        double max_r = fabs(r[1]);
        for (int i = 2; i <= n; i++) {
            if (fabs(r[i]) > max_r) {
                max_r = fabs(r[i]);
                k = i;
            }
        }

        if (max_r > EPS) stop = 1;

        if (stop) {
            double d = r[k];
            x[k] += d;
            for (int i = 1; i <= n; i++) {
                r[i] -= a[i][k] * d;
            }
            lan_lap++;
        }
        if (lan_lap > 500) break;
    } while (stop);

    printf("\n+------+");
    for(int i=1; i<=n; i++) printf("------------+");
    for(int i=1; i<=n; i++) printf("------------+");
}

// NỘI SUY LAGRANGE
void NoiSuyLagrange(double x[], double y[], int n, double X) {
    double result = 0;
    for (int i = 0; i < n; i++) {
        double pi = 1;
        for (int j = 0; j < n; j++) {
            if (j != i) pi *= ((X - x[j]) / (x[i] - x[j]));
        }
        result += y[i] * pi;
    }
    printf("\n=== KET QUA NOI SUY LAGRANGE ===\n");
    printf("Gia tri f(%.3f) = %.6f\n", X, result);
}

// BÌNH PHƯƠNG BÉ NHẤT
void HoiQuyMenu(double x[], double y[], int n) {
    int lc;
    printf("\n======= HOI QUY =======\n");
    printf("1. Tuyến tính: y = a + bx\n");
    printf("2. Bậc hai:    y = a + bx + cx^2\n");
    printf("3. Hàm mũ:     y = a * e^(bx)\n");
    printf("4. Lũy thừa:   y = a * x^b\n");
    printf("Nhập lựa chọn loại hình: ");
    scanf("%d", &lc);

    in_bang_hoi_quy(x, y, n, (lc == 2));

    if (lc == 1) {
        double tx = tinh_mu(x,1,y,0,n), ty = tinh_mu(x,0,y,1,n);
        double tx2 = tinh_mu(x,2,y,0,n), txy = tinh_mu(x,1,y,1,n);
        double b = (txy - (ty * tx) / n) / (tx2 - (tx * tx) / n);
        double a = (ty - b * tx) / n;
        printf("\n=> Phuong trinh tim duoc: y = %f + %f * x\n", a, b);
    }
    else if (lc == 2) {
        double tx = tinh_mu(x,1,y,0,n), tx2 = tinh_mu(x,2,y,0,n);
        double tx3 = tinh_mu(x,3,y,0,n), tx4 = tinh_mu(x,4,y,0,n);
        double ty = tinh_mu(x,0,y,1,n), txy = tinh_mu(x,1,y,1,n), tx2y = tinh_mu(x,2,y,1,n);

        double D = n*(tx2*tx4 - tx3*tx3) - tx*(tx*tx4 - tx2*tx3) + tx2*(tx*tx3 - tx2*tx2);
        double Da = ty*(tx2*tx4 - tx3*tx3) - tx*(txy*tx4 - tx3*tx2y) + tx2*(txy*tx3 - tx2*tx2y);
        double Db = n*(txy*tx4 - tx3*tx2y) - ty*(tx*tx4 - tx2*tx3) + tx2*(tx*tx2y - txy*tx2);
        double Dc = n*(tx2*tx2y - txy*tx3) - tx*(tx*tx2y - ty*tx3) + ty*(tx*txy - tx2*tx2);

        if (fabs(D) < 1e-9) return;
        printf("\n=> Phuong trinh tim duoc: y = %f + %f * x + %f * x^2\n", Da/D, Db/D, Dc/D);
    }
    else if (lc == 3 || lc == 4) {
        double ln_x[MAX], ln_y[MAX];
        for (int i = 0; i < n; i++) {
            if (y[i] <= 0 || (lc == 4 && x[i] <= 0)) return;

            ln_x[i] = log(x[i]);
            ln_y[i] = log(y[i]);
        }
        if (lc == 3) { // Hàm mũ
            double tx = tinh_mu(x,1,ln_y,0,n), t_lny = tinh_mu(x,0,ln_y,1,n);
            double tx2 = tinh_mu(x,2,ln_y,0,n), tx_lny = tinh_mu(x,1,ln_y,1,n);
            double b = (tx_lny - (t_lny * tx) / n) / (tx2 - (tx * tx) / n);
            double a = exp((t_lny - b * tx) / n);
            printf("\n=> Phuong trinh tim duoc: y = %f * e^(%f * x)\n", a, b);
        } else { // Hàm lũy thừa
            double t_lnx = tinh_mu(ln_x,1,ln_y,0,n), t_lny = tinh_mu(ln_x,0,ln_y,1,n);
            double t_lnx2 = tinh_mu(ln_x,2,ln_y,0,n), t_lnxlny = tinh_mu(ln_x,1,ln_y,1,n);
            double b = (t_lnxlny - (t_lny * t_lnx) / n) / (t_lnx2 - (t_lnx * t_lnx) / n);
            double a = exp((t_lny - b * t_lnx) / n);
            printf("\n=> Phuong trinh tim duoc: y = %f * x^(%f)\n", a, b);
        }
    }
}

int doc_file_he_pt(const char *filename, double a[MAX][MAX]) {
    FILE *f = fopen(filename, "r");
    if (!f) { printf("Khong the mo file %s\n", filename); return -1; }
    int n;
    if (fscanf(f, "%d", &n) != 1) { fclose(f); return -1; }
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n + 1; j++) {
            if (fscanf(f, "%lf", &a[i][j]) != 1) { fclose(f); return -1; }
        }
    }
    fclose(f);
    return n;
}

// ĐỌC FILE NỘI SUY
int doc_file_mocsuy(const char *filename, double x[], double y[], double *X) {
    FILE *f = fopen(filename, "r");
    if (!f) { printf("Khong the mo file %s\n", filename); return -1; }
    int n;
    if (fscanf(f, "%d", &n) != 1) { fclose(f); return -1; }
    for (int i = 0; i < n; i++) fscanf(f, "%lf", &x[i]);
    for (int i = 0; i < n; i++) fscanf(f, "%lf", &y[i]);
    if (X != NULL) {
        if (fscanf(f, "%lf", X) != 1) *X = 0;
    }
    fclose(f);
    return n;
}

int main(){
    int dieukhien = -1;
    double A[MAX][MAX], x[MAX], y[MAX], X_target;
    int n;
    char ten_file[100];

    while (dieukhien != 0){
        printf("\n========MENU CHUONG TRINH==========\n");
        printf("1. Giai he PT tuyen tinh bang Cramer\n");
        printf("2. Giai he PT tuyen tinh bang Gauss\n");
        printf("3. Giai he PT tuyen tinh bang Gauss-Seidel\n");
        printf("4. Giai he PT tuyen tinh bang Giam du\n");
        printf("5. Noi suy Lagrange\n");
        printf("6. Noi suy Binh phuong be nhat\n");
        printf("0. Thoat chuong trinh\n");
        printf("========================================\n");
        printf("Nhap lua chon: ");
        scanf("%d", &dieukhien);

        if (dieukhien >= 1 && dieukhien <= 6){
            printf("Nhap ten file: ");
            scanf("%s", ten_file);
        }

        switch (dieukhien){
            case 1:
                n = doc_file_he_pt(ten_file, A);
                if (n > 0) Cramer(A, n);
                break;

            case 2:
                n = doc_file_he_pt(ten_file, A);
                if (n > 0) Gauss(A, n);
                break;

            case 3:
                n = doc_file_he_pt(ten_file, A);
                if (n > 0) GaussSeidel(A, n);
                break;

            case 4:
                n = doc_file_he_pt(ten_file, A);
                if (n > 0) GiamDu(A, n);
                break;

            case 5:
                n = doc_file_mocsuy(ten_file, x, y, &X_target);
                if (n > 0) NoiSuyLagrange(x, y, n, X_target);
                break;

            case 6:
                n = doc_file_mocsuy(ten_file, x, y, NULL);
                if (n > 0) HoiQuyMenu(x, y, n);
                break;
            case 0:
                printf("\nTam biet!\n");
                break;

            default:
                printf("\nLua chon khong hop le, vui len chon lai!\n");
        }
    }
    return 0;
}
