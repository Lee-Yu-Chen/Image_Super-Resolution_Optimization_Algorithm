
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
// 
// 
// 15 December 2023
// 
// - Brief discussion of the project plan in study break
// - Update on current progress : in the middle of implementing the mathematical concept of the ML-like algo
// - Pro forma, project review meeting and meeting minutes
// 
// 
// 
// 
// 
// 5 Jan 2023
// 
// Agenda : 
// - Explanation of method of implementation of finding coefficients and image reconstruction
// - Demonstration of result, i.e. image reconstruction
// - Review of the quality of reconstructed image and values of coefficients
// 
// Next things to do : 
// - Use a more advance method of implementing the upsampling and finding coeffcient : 
//      - instead of upsampling by 2, umsample by 10
//      - each of the intermediate pixels take references from 9 surrounding pixels from low-res image
// 
// - Start MISR
// - Finding the misalignment between 2 images by fitting a 2-dimensional parabolic curve to each image, 
// - Use what called "Design of Experiments" method to find the coefficients of the parabolic curves
// - then find the distance between the peaks of the two 2-D parabolic curves, and this is where the 2 images correlate the most,
// - then x and y offset can be obtained, then the 2 images can be aligned  
// 
// 
// Notes : 
// - thesis need to include Math
// - Pro forma
// - Image rotation is too complicated for this, maybe suitable for a PhD
// - Correlation should be non-integer
// - It is impossible to get image better than the reference image (the original hi-res image) unless we're using MISR
// - Obtaining significant improvement from the low-res image is already very good 
// - Just keep going and we'll see how far we can go forward, 
// - until one day we'll need to stop and just converge what we've got and write them into the thesis
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 7 Feb 2024
// 
// - check through the R Matrix genaration, suspect bug in the R matrix generation bc shouldn't suppose to be 0 determinance
// - suspect dont have enough summation 
// 
// 
// - try using a small dummy image to check
// 
// potential problems:
// - got bug in R matrix generatrion
// - overflowing
// - mis-define of R matrix  mis-definition)
//
//
//
//
//
//  
// 
// 
// 
// What I did : 
// - used a cropped image
// - change 255 to 65535
// - change x0x0-x8x8 long long int to unsigned long long int 
// 
// - change x8x7 back to long long int : : overflow and debug
// 
// - check Arr[C] : correct
// - check x0-x8 indexing : correct
// 
// - debug 49601*49601 overflowing problem
// 
// - change x0-x8 from int to long long int
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
//  21 Feb 2024
//
// Comment : 
// - The result wasn't very good, but it is expected and generally acceptable, as only 1 image is used for reconstrcution
// - Also, the downsampling-upsampling scale was 5, means that it is taking 1/25 of the information from the low-res image, and in this case certain degree of blurriness is expected
// - Also, due to the nature of the method used to reconstruct the image, artifacts that exist in the reconstructed image is also expected
// 
//
// Next :
// -correlation and interpolation
// - 
//
//
// How to get test images : 
// - use tripod or any stable camera, dont move the camera and shoot multiple photos
// - use a scene where there are no repetitive patterns, not a large blank surface, and better if there are some detailed patterns like words to be correlated and can be recognised by human eyes
// - use a digital way : down sample image with different indexing
// 
// 
// 
// Thing :
// - I need test images/videos for MISR
// - Pro forma
// 
// 
// 
//
// 28 Feb 2024
// 
// - choose a smaller vicinity data for the curve fitting, otherwise it would not be a good fit
// - a small vicinity of 12*12 or 16*16 around the highest data point would be good
// - 
// - plot the sampled data and curve on the same plane to check if it is a good fit
// - differentiate with respect to x and y to obtain the peak of the bipolar curve (dz/dy and dz/dx)
// - MATLAB command 
// x = -3:0.1:3;
// y = -3:0.1:3;
// [xx,yy] = meshgrid(x,y)
// 
// 
// 
// 
// 
// 
// 6 March 2024
// 
// - Demonstration of the curve fitted to a small vicinity of size 21, 7
// - Discussion about the reasons of not getting a good fit
// - Attempting curve fitting with a smaller vicinity of size 3, debug and found out the problem causing a weird curve, 
// - turn out the curve fitting result with vicinity size of 3 is satisfactory
// 
// Next:
// 
// - can prepare to move on to the next stage of MISR i.e. correlation + advanced ML-like algo
// - the exact methodology and math of the next stage of MISR have not been sorted out yet
// - try the correlation with 2 different but close-to-each-other photos and verify it is working i.e. can tell how much the image is shifted relatively
// 
// 
// 
// 
// 
// 
// 
// 13 March 2024 (Meeting minute)
// 
// - complete the correlation, finding offset, interpolation method first
// - complete the objectives as much as possible and make it looks wrapped up
// - (because Dr Mumtaj and the other 2 academics that are going to assess me would probably assess me based on the project objectives and deliverables)
// - do thesis in parallel and see if have time for MISR
// 
// - use cropped images of pseudo-identical taken photos as test images
// - try AOI_size = 5
// 
// - maximization / minimization => optimization
// - DOE (design of experiment)
// 
// 
// 
// 
// 
// 
// 
// 19-20 March 2024 (log)
// 
// - determined the x and y position of the peak of the fitted 2D-curve to the correlation plot (the correlation array)
// - from the x and y position of the peak, determined the x and y offset (the alignment) of the 2 images
// - generate an image from the 2 images by taking the average value of each pixel in the images to verify the value of the offsets found 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 




/*


#include <Magick++.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <Eigen/Dense>

using namespace std;
using namespace Magick;
//using Eigen::MatrixXi;


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



    const int s = 5;      // scaling factor



    Image image;
    Image image2;
    //image.read("C:/Users/ACER/Downloads/Picture/Batu Ferringhi.jpg");

    //image.read("C:/Users/ACER/Desktop/C/C++/Y4 FYP/Test Images/test.png");
    //Image image2 = image;

    image.read("C:/Users/ACER/Desktop/C/C++/Y4 FYP/Test Images/IMG_20240305_172658.jpg");
    image2.read("C:/Users/ACER/Desktop/C/C++/Y4 FYP/Test Images/IMG_20240305_172659.jpg");



    int dimen = 500;            // dimensions of the cropped images
    int x_crop_offs = 850;        // x and y offsets for cropping the images
    int y_crop_offs = 1400;      // dimen = 480, x_crop_offs = 20, y_crop_offs = 1100 for test image (testing MISR with a single image test.png)

    int x_offs = -14;
    int y_offs = -19;            // manual x and y offsets


    image.crop(Geometry(dimen, dimen, x_crop_offs, y_crop_offs));
    image2.crop(Geometry(dimen, dimen, x_crop_offs + x_offs, y_crop_offs + y_offs));

    //image2.roll(120, 200);




    size_t w = image.size().width();
    size_t h = image.size().height();

    size_t w2 = image2.size().width();
    size_t h2 = image2.size().height();



    //image.sample(Geometry(w/s, h/s));
    image.write("C:/Users/ACER/Desktop/asd.png");
    //image.read("C:/Users/ACER/Desktop/sampled.png");

    image2.write("C:/Users/ACER/Desktop/asd2.png");



    /*
    Image image;
    image.read("C:/Users/ACER/Downloads/Picture/Aqua.jpg");
    size_t width = image.size().width()/2;
    size_t height = image.size().height()/2;
    image.sample(Geometry(width,height));
    *




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
    Image im;
    im.read("C:/Users/ACER/Desktop/asd.png");

    size_t width = im.size().width();
    size_t height = im.size().height();

    //int m[width][height]          this only valid in C, but not in C++

    //int* img1 = (int*)calloc(width * height * 3, sizeof(int));
    //int* img2 = (int*)calloc(width * height * 3, sizeof(int));

    int* img1 = new int[width * height * 3];
    int* img2 = new int[width * height * 3];


    Pixels pixel(image);
    Pixels pixel2(image2);

    Quantum* pix = pixel.get(0, 0, width, height);
    Quantum* pix2 = pixel2.get(0, 0, w2, h2);


    // moving the range of pixel values from 0~65535 to -32768~32767
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




    double* X = new double[(width * 2 - 1) * (height * 2 - 1)];      // a fitted 2D curve, supposed to be very smooth
    double* Y = new double[(width * 2 - 1) * (height * 2 - 1)];      // actual value sampled in the mesh() in MATLAB, aka correlation array


    I = 0;

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
                            getPixelRGB(image2, n + width - 1 - i, m + height - 1 - j, R2, G2, B2);  // ok

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

            ///////////////////////
            sum = Rsum + Gsum + Bsum;                   // this line is doubtful

            corr = (double)sum / (double)Nsum;
            //corr = (double)Rsum / (double)Nsum;


            //append("C:/Users/ACER/Desktop/C/C++/Y4 FYP/test folder/value.txt", to_string(corr) + " ");      // this line also takes quite long time as open file muitlple times
            write << to_string(corr) + " ";

            //cout << i << "  ,  " << j << "\t\t" << Nsum << '\t' << sum << endl;

            *(Y + I) = corr;
            I++;
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




    long long int highest = -999999999999999;
    int x_posi;
    int y_posi;

    int arr_width = 2 * width - 1;
    int arr_height = 2 * height - 1;

    //for (int i = 0; i <= 2 * height - 2; i++) {
        //for (int j = 0; j <= 2 * width - 2; j++) {

    int width_crop = (w + w2) / 10;;      // area to be cropped out in the correlation array, basically the edges
    int height_crop = (h + h2) / 10;

    for (int i = height_crop; i <= 2 * height - height_crop - 2; i++) {
        for (int j = width_crop; j <= 2 * width - width_crop - 2; j++) {

            I = i * arr_width + j;

            if (*(Y + I) > highest) {
                highest = *(Y + I);
                x_posi = j + 1;
                y_posi = i + 1;
            }

        }
    }


    cout << highest << '\t' << x_posi << '\t' << y_posi << endl;






    // extracting the array of interest (AOI) from the correlation array
    // highest sample located at 333,380
    // range of AOI for size 21 is 323-343, 370-390
    // range of AOI for size 7 is 330-336, 377-383

    int AOI_size = 5;                                           // should be odd number
    double* AOI = new double[AOI_size * AOI_size];             // array of interest

    //int start_row = 333 - (AOI_size - 1) / 2 - 1;
    //int start_column = 380 - (AOI_size - 1) / 2 - 1;
    int start_row = y_posi - (AOI_size - 1) / 2 - 1;
    int start_column = x_posi - (AOI_size - 1) / 2 - 1;

    int rows = (height * 2) - 1;
    int columns = (width * 2) - 1;

    ofstream write1("C:/Users/ACER/Desktop/C/C++/Y4 FYP/test folder/Array of Interest.txt");
    write1 << endl;


    // collecting samples into the AOI array
    for (int i = 0; i <= AOI_size - 1; i++) {                 // row
        for (int j = 0; j <= AOI_size - 1; j++) {             // column

            *(AOI + i * AOI_size + j) = *(Y + (start_row + i) * columns + (start_column + j));             // this line verified through MATLAB
        }
    }


    // writting AOI sample data from array into txt file
    for (int i = 0; i <= AOI_size - 1; i++) {
        for (int j = 0; j <= AOI_size - 1; j++) {

            cout << *(AOI + i * AOI_size + j) << "  ";
            write1 << to_string(*(AOI + i * AOI_size + j)) + " ";
        }
        cout << endl;
        write1 << ";\n";
    }

    write1.close();

    // 21
    // z = -834690.*xx.*xx + 16557900.*xx - 121994.*yy.*yy + 2303980.*yy + 13590.3.*xx.*yy + 512008000;
    // 
    // 7
    // z = -6545100. * xx.*xx + 39175600. * xx - 4364430. * yy.*yy + 26091600. * yy + 31673.5.*xx.*yy + 571094000;



    // z = coeff(1)*xx.*xx + coeff(2)*xx + coeff(3)*yy.*yy + coeff(4)*yy + coeff(5)*xx.*yy + coeff(6);
    // z = coeff13(1)*xx.*xx + coeff13(2)*xx + coeff13(3)*yy.*yy + coeff13(4)*yy + coeff13(5)*xx.*yy + coeff13(6);
    // mesh(xx,yy,z)


    // z = coeff(1)*yy.*yy + coeff(2)*yy + coeff(3)*xx.*xx + coeff(4)*xx + coeff(5)*xx.*yy + coeff(6);


























    //////////////////////////////////////
    // finding coefficient start here



    Eigen::Matrix<long long int, 6, 6> R;       // R matrix of the curve fitting optimisation

    Eigen::Matrix<double, 6, 6> R_double;

    Eigen::Matrix<double, 6, 1> P, curve_coeff;           // P matrix and the coefficient matrix (A matrix) of the 2D curve fitting optimisation


    R << 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0;

    P << 0, 0, 0, 0, 0, 0;

    I = 0;

    double y;

    cout << "\n\nstarting fitting curve and finding coefficients...\n\n";
    //for (int x0 = 0; x0 <= (height - 1) * 2; x0++) {            // let x0 = row
      // for (int x1 = 0; x1 <= (width - 1) * 2; x1++) {         // let x1 = column


    //for (int x0 = 0; x0 <= AOI_size - 1; x0++) {             // let x0 = row
        //for (int x1 = 0; x1 <= AOI_size - 1; x1++) {         // let x1 = column


    for (int x0 = 1; x0 <= AOI_size; x0++) {             // let x0 = row
        for (int x1 = 1; x1 <= AOI_size; x1++) {         // let x1 = column
            //for (int x0 = x_posi-1; x0 <= x_posi+1; x0++) {             // let x0 = row
                //for (int x1 = y_posi-1; x1 <= y_posi+1; x1++) {         // let x1 = column

                    //r, c
            R(0, 0) += x0 * x0 * x0 * x0;
            R(1, 0) += x0 * x0 * x0;
            R(2, 0) += x0 * x0 * x1 * x1;
            R(3, 0) += x0 * x0 * x1;
            R(4, 0) += x0 * x0 * x0 * x1;
            R(5, 0) += x0 * x0;


            R(0, 1) += x0 * x0 * x0;
            R(1, 1) += x0 * x0;
            R(2, 1) += x0 * x1 * x1;
            R(3, 1) += x0 * x1;
            R(4, 1) += x0 * x0 * x1;
            R(5, 1) += x0;


            R(0, 2) += x0 * x0 * x1 * x1;
            R(1, 2) += x0 * x1 * x1;
            R(2, 2) += x1 * x1 * x1 * x1;
            R(3, 2) += x1 * x1 * x1;
            R(4, 2) += x0 * x1 * x1 * x1;
            R(5, 2) += x1 * x1;


            R(0, 3) += x0 * x0 * x1;
            R(1, 3) += x0 * x1;
            R(2, 3) += x1 * x1 * x1;
            R(3, 3) += x1 * x1;
            R(4, 3) += x0 * x1 * x1;
            R(5, 3) += x1;


            R(0, 4) += x0 * x0 * x0 * x1;
            R(1, 4) += x0 * x0 * x1;
            R(2, 4) += x0 * x1 * x1 * x1;
            R(3, 4) += x0 * x1 * x1;
            R(4, 4) += x0 * x0 * x1 * x1;
            R(5, 4) += x0 * x1;


            R(0, 5) += x0 * x0;
            R(1, 5) += x0;
            R(2, 5) += x1 * x1;
            R(3, 5) += x1;
            R(4, 5) += x0 * x1;
            R(5, 5) += 1;




            //y = *(Y + I);
            y = *(AOI + I);
            I++;

            P(0) += x0 * x0 * y;
            P(1) += x0 * y;
            P(2) += x1 * x1 * y;
            P(3) += x1 * y;
            P(4) += x0 * x1 * y;
            P(5) += y;

            //cout << R(0, 0) << endl;
        }
    }


    R_double << R(0, 0), R(0, 1), R(0, 2), R(0, 3), R(0, 4), R(0, 5),
        R(1, 0), R(1, 1), R(1, 2), R(1, 3), R(1, 4), R(1, 5),
        R(2, 0), R(2, 1), R(2, 2), R(2, 3), R(2, 4), R(2, 5),
        R(3, 0), R(3, 1), R(3, 2), R(3, 3), R(3, 4), R(3, 5),
        R(4, 0), R(4, 1), R(4, 2), R(4, 3), R(4, 4), R(4, 5),
        R(5, 0), R(5, 1), R(5, 2), R(5, 3), R(5, 4), R(5, 5);


    cout << "\n\n" << R << "\n\n";

    cout << "\n\n" << R_double << "\n\n\n\n";


    cout << P << "\n\n\n\n";


    curve_coeff = R_double.inverse() * P;

    cout << curve_coeff << "\n\n";

    double a = curve_coeff(0);
    double b = curve_coeff(1);
    double c = curve_coeff(2);
    double d = curve_coeff(3);
    double e = curve_coeff(4);

    double x_ = (2 * b * c - d * e) / (e * e - 4 * a * c);
    double y_ = 2 * a * (d * e - 2 * b * c) / (e * (e * e - 4 * a * c)) - b / e;


    cout << "\nx_  " << x_ << "\ty_  " << y_ << endl;

    int minus = AOI_size / 2 + 1;

    cout << minus << endl << endl;


    double x_offset = x_ - minus + x_posi - dimen;
    double y_offset = y_ - minus + y_posi - dimen;

    cout << "x_offset  " << x_offset << endl;
    cout << "y_offset  " << y_offset << endl;











    ///////////////////////////////
    // combining two images to verify the x and y offsets found in the correlation method



    // generate a image by combining asd.png and asd2.png 
    Image combined(Geometry(width, height), Color("#000000"));
    Pixels p(combined);
    Quantum* q = p.get(0, 0, width, height);

    cout << "\n\nCombining ...\n\n";

    I = 0;
    for (int i = 0; i <= height - 1; i++) {
        for (int j = 0; j <= width - 1; j++) {

            *(q + I) = (*(pix + I) + *(pix2 + I)) / 2;
            *(q + I + 1) = (*(pix + I + 1) + *(pix2 + I + 1)) / 2;
            *(q + I + 2) = (*(pix + I + 2) + *(pix2 + I + 2)) / 2;
            I += 3;

            //cout << I << endl;

        }
    }

    p.sync();
    combined.write("C:/Users/ACER/Desktop/combined.png");











































    //free(img1);
    //free(img2);
    delete img1;
    delete img2;

    delete Y;
    delete X;



    return 0;
}




*/



