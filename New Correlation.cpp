
// things to add : 
// - think of a method to do preprocessing of the normalisation - move the minus outside of the for loop
// - find out the reason and solution to the low CPU utilisation
// - change the logic flow -> for for for for if else, if possible
// - debugging for unsuccessful colours
// - getPixelRGB() can also be optimised : pixel.get() a large area and check pixel one by one, but not declare Quantum* pixel and pixel.get() every time
// - to make the program fast, do as much pre-processing as possible


// progress review meeting & pro forma
// meeting in study break
// **correlation as one milestone

// changes made:
// - develop a new method to access pixels
// - remove usage of function
// - data type long int -> long long int
// - float -> double
// 
// data overflowing problem




// **make it take 2 images with different sizes
// **make it more automatic

// 8 Dec 2023
// Demonstration of Cost Function.cpp
// Discussion of problem overcame while developing the program ( int, long int overflow problem)
// inspection of memory allocation for different data type as verification of the problem identified
// discussion of math theory behind the weighted sum/ML-like algorithm

// 11 Dec 2023
// 
// Demonstration of New Correlation.cpp
// dicussion of improvements to the program
//  - can be used for images with different sizes
//  - increase the level of automation in the program ( after collecting data in file, use a MATLAB script to plot it/use a C++ code to find the peaks)
// 
// the direction of the upcoming tasks : usage of C++ libraries to implement the mathematics of the ML-like algo
// the direction of the future implementation
// 
// 
// meeting in the study break
// progress review meeting
//
//

/*



#include <Magick++.h>
#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;
using namespace Magick;


void append(string txtfile, string text) {

    // reading the exising contents in the file
    fstream read(txtfile);

    // store the file contents in a variable
    string storage, buffer;
    while (getline(read, buffer)) {
        storage += buffer;
        storage += '\n';
    }

    // erase the last '\n'
    storage.erase(storage.length() - 1, 1);

    // add the new text at the end of the file
    storage += text;
    read.close();

    fstream write(txtfile);

    // write the existing contents + the added text to the file
    write << storage;
    write.close();
};


//void getPixelRGB(Image img, int x, int y, int* R, int* G, int* B) {
void getPixelRGB(Image img, int x, int y, Quantum* R, Quantum* G, Quantum* B) {
    Pixels pixel(img);

    Quantum* pix = pixel.get(x, y, 1, 1);

    *R = *pix;
    *G = *(pix + 1);
    *B = *(pix + 2);
};







int main(int argc, char** argv)
{
    InitializeMagick(*argv);

    /*
    // create a 16x16 white image with a red dot
    Image image(Geometry(16, 16), Color("#000000"));         // #7df9ff
    Pixels p(image);
    Quantum* q = p.get(3, 3, 1, 1);

    *q = 65535;
    *(q + 1) = 0;
    *(q + 2) = 0;

    p.sync();
    image.write("C:/Users/ACER/Desktop/asd.png");
    *



    Image image;
    //image.read("C:/Users/ACER/Downloads/Picture/Batu Ferringhi.jpg");
    //image.write("C:/Users/ACER/Desktop/Original.png");
    image.read("C:/Users/ACER/Desktop/sampled.png");



    /*
    Image image;
    image.read("C:/Users/ACER/Downloads/Picture/Aqua.jpg");
    size_t width = image.size().width()/2;
    size_t height = image.size().height()/2;
    image.sample(Geometry(width,height));
    *








    Image rolled = image;
    rolled.roll(120, 200);
    rolled.write("C:/Users/ACER/Desktop/rolled.png");


    ofstream write("C:/Users/ACER/Desktop/C/C++/Y4 FYP/test folder/value.txt");
    write << endl;
    //write.close();









    /*
    // test code for getPixelRGB()
    Quantum r1, r2, g1, g2, b1, b2;
    Quantum* R1 = &r1;
    Quantum* R2 = &r2;
    Quantum* G1 = &g1;
    Quantum* G2 = &g2;
    Quantum* B1 = &b1;
    Quantum* B2 = &b2;

    int r1, r2, g1, g2, b1, b2;
    int* R1 = &r1;
    int* R2 = &r2;
    int* G1 = &g1;
    int* G2 = &g2;
    int* B1 = &b1;
    int* B2 = &b2;


    size_t wid = image.size().width();
    size_t hei = image.size().height();
    for (int i = 0; i < wid; i++) {
        for (int j = 0; j < hei; j++) {

            getPixelRGB(image, i, j, R1, G1, B1);

            cout << i << '\t' << j << '\t' << *R1 << '\t' << *G1 << '\t' << *B1 << endl;

        }
    }
    *





    // need to add '\n' at the front of the added text
    //append("C:/Users/ACER/Desktop/C/C++/Y4 FYP/test folder/value.txt", " Helo");

    // image dimensions
    // assumming two images must have same dimensions

    size_t width = image.size().width();
    size_t height = image.size().height();

    //int m[width][height]          this only valid in C, but not in C++

    //int* img1 = (int*)calloc(width * height * 3, sizeof(int));
    //int* img2 = (int*)calloc(width * height * 3, sizeof(int));

    int* img1 = new int[width * height * 3];
    int* img2 = new int[width * height * 3];


    Pixels pixel(image);
    Pixels pixel2(rolled);

    Quantum* pix = pixel.get(0, 0, width, height);
    Quantum* pix2 = pixel2.get(0, 0, width, height);

    int I = -1;
    for (int i = 0; i <= height - 1; i++) {
        for (int j = 0; j <= width - 1; j++) {
            for (int k = 0; k <= 2; k++) {
                I++;
                *(img1 + I) = (*(pix + I)) - (65536 / 2);
                *(img2 + I) = (*(pix2 + I)) - (65536 / 2);
            }
        }
        //cout << i << endl;
    }
    cout << "done allocating" << endl;

    /*
    // test code to examine the allocated memory array
    Image test(Geometry(width, height), Color("#000000"));
    Pixels p(test);
    Quantum* q = p.get(0, 0, width, height);

    int J=-3;
    for (int i = 0; i <= height - 1; i+=1) {
        for (int j = 0; j <= width - 1; j+=2) {
            J = ((i * width) + j) * 3;
            //J += 3;
            *(q + J) = *(img1 + J);
            *(q + J + 1) = *(img1 + J + 1);
            *(q + J + 2) = *(img1 + J + 2);
        }
        //cout << i << endl;
    }
    cout << "done re-creating" << endl;

    p.sync();
    test.write("C:/Users/ACER/Desktop/test.png");
    *






    long long int Rsum, Gsum, Bsum;      // change to long long int?
    double corr;
    long long int Nsum, sum;
    Quantum r1, r2, g1, g2, b1, b2;
    Quantum* R1 = &r1;
    Quantum* R2 = &r2;
    Quantum* G1 = &g1;
    Quantum* G2 = &g2;
    Quantum* B1 = &b1;
    Quantum* B2 = &b2;

    time_t t;
    char tim[26] = {};



    int J1, J2;
    int x1, x2, y1, y2;

    int wm1 = width - 1;
    int hm1 = height - 1;





    // for for for for if else


    //// using img1 coordinate as reference
    //// j : y position of img2 with respect to img1
    //// i : x position of img2 with respect to img2


    // definition
    // i : how much img1 is higher than img2, 0 -> img1 lower, 99 -> img1 higher
    // j : how much img1 is lefter than img2, 0 -> img1 right, 99 -> img1 left
    // i,j are alignments
    // m,n are pixel coordinate
    for (int j = 0; j <= (height - 1) * 2; j++) {

        //append("C:/Users/ACER/Desktop/C/C++/Y4 FYP/test folder/value.txt", "\n\n");

        for (int i = 0; i <= (width - 1) * 2; i++) {

            //cout << "start computing" << endl;

            // for every alignment

            sum = 0;
            Nsum = 0;
            Rsum = 0;
            Gsum = 0;
            Bsum = 0;


            if (j <= height - 1) {

                for (int m = 0; m <= j; m++) {

                    if (i <= width - 1) {         // i from 0 to 99

                        for (int n = 0; n <= i; n++) {                  // from short to long, start with 0

                            // ascending m, ascending n
                            /*
                            getPixelRGB(image, n, m, R1, G1, B1);
                            getPixelRGB(rolled, n + width - 1 - i, m + height - 1 - j, R2, G2, B2);  // ok

                            Rsum += (*R1 - (65536 / 2)) * (*R2 - (65535 / 2));          // move outside of the for loop
                            Gsum += (*G1 - (65536 / 2)) * (*G2 - (65535 / 2));
                            /Bsum += (*B1 - (65536 / 2)) * (*B2 - (65535 / 2));
                            */

                            /*
                            x1 = n;
                            y1 = m;
                            x2 = n + width - 1 - i;
                            y2 = m + height - 1 - j;

                            J1 = ((m * width) + n) * 3;
                            J2 = (((m + hm1 - j) * width) + (n + wm1 - i)) * 3;
                            *

                            Rsum += (*(img1 + (((m * width) + n) * 3))) * (*(img2 + ((((m + hm1 - j) * width) + (n + wm1 - i)) * 3)));
                            Gsum += (*(img1 + (((m * width) + n) * 3) + 1)) * (*(img2 + ((((m + hm1 - j) * width) + (n + wm1 - i)) * 3) + 1));
                            Bsum += (*(img1 + (((m * width) + n) * 3) + 2)) * (*(img2 + ((((m + hm1 - j) * width) + (n + wm1 - i)) * 3) + 2));

                            Nsum++;

                        }
                    }

                    else if (i > width - 1) {     // i from 100 to 198

                        for (int n = width - 1; n >= i - width + 1; n--) {    // from long to short, start with 99

                            // ascending m, descending n
                            Rsum += (*(img1 + (((m * width) + n) * 3))) * (*(img2 + ((((m + hm1 - j) * width) + (n + wm1 - i)) * 3)));
                            Gsum += (*(img1 + (((m * width) + n) * 3) + 1)) * (*(img2 + ((((m + hm1 - j) * width) + (n + wm1 - i)) * 3) + 1));
                            Bsum += (*(img1 + (((m * width) + n) * 3) + 2)) * (*(img2 + ((((m + hm1 - j) * width) + (n + wm1 - i)) * 3) + 2));

                            Nsum++;
                        }
                    }
                }
            }
            else if (j > height - 1) {

                for (int m = height - 1; m >= j - height + 1; m--) {

                    if (i <= width - 1) {

                        for (int n = 0; n <= i; n++) {

                            // descending m, ascending n
                            Rsum += (*(img1 + (((m * width) + n) * 3))) * (*(img2 + ((((m + hm1 - j) * width) + (n + wm1 - i)) * 3)));
                            Gsum += (*(img1 + (((m * width) + n) * 3) + 1)) * (*(img2 + ((((m + hm1 - j) * width) + (n + wm1 - i)) * 3) + 1));
                            Bsum += (*(img1 + (((m * width) + n) * 3) + 2)) * (*(img2 + ((((m + hm1 - j) * width) + (n + wm1 - i)) * 3) + 2));

                            Nsum++;
                        }
                    }
                    else if (i > width - 1) {

                        for (int n = width - 1; n >= i - width + 1; n--) {

                            // descending m, descending n
                            Rsum += (*(img1 + (((m * width) + n) * 3))) * (*(img2 + ((((m + hm1 - j) * width) + (n + wm1 - i)) * 3)));
                            Gsum += (*(img1 + (((m * width) + n) * 3) + 1)) * (*(img2 + ((((m + hm1 - j) * width) + (n + wm1 - i)) * 3) + 1));
                            Bsum += (*(img1 + (((m * width) + n) * 3) + 2)) * (*(img2 + ((((m + hm1 - j) * width) + (n + wm1 - i)) * 3) + 2));

                            Nsum++;
                        }
                    }
                }
            }


            sum = Rsum + Gsum + Bsum;

            corr = (double)sum / (double)Nsum;
            //corr = (double)Rsum / (double)Nsum;


            //append("C:/Users/ACER/Desktop/C/C++/Y4 FYP/test folder/value.txt", to_string(corr) + " ");      // this line also takes quite long time as open file muitlple times
            write << to_string(corr) + " ";

            //cout << i << "  ,  " << j << "\t\t" << Nsum << '\t' << sum << endl;
        }
        //append("C:/Users/ACER/Desktop/C/C++/Y4 FYP/test folder/value.txt", ";");
        write << ";\n";

        t = time(0);
        ctime_s(tim, 26, &t);
        cout << j << '\t' << tim << endl;
    }

    // started 3.44am 30 Nov 2023
    // only reach row 18 (alignment) by 4.06am
    // using aroung 25% of CPU usage



     //getPixelRGB(test, i, j, R, G, B);



    write.close();


    //free(img1);
    //free(img2);
    delete img1;
    delete img2;



    /*
    ofstream write("C:/Users/ACER/Desktop/C/C++/Y4 FYP/test folder/test.txt");

    write << "asdasd";
    write << '\n';
    //write << "1234567";
    //write << '\t';
    write << "i love you";

    write.close();
    *



    return 0;
}




*/




