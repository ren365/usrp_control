// set up our masks, defining the pin numbers
#define AMP_GPIO_MASK   (1 << 6) | (1 << 10)
#define MAN_GPIO_MASK   (1 << 4)
#define ATR_MASKS       (AMP_GPIO_MASK | MAN_GPIO_MASK)
// set up our values for ATR control: 1 for ATR, 0 for manual
#define ATR_CONTROL     (AMP_GPIO_MASK)
// set up the GPIO directions: 1 for output, 0 for input
#define GPIO_DDR        (AMP_GPIO_MASK | MAN_GPIO_MASK)
// assume an existing multi_usrp object, called "usrp"
// now, let's do the basic ATR setup
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/safe_main.hpp>
#include <iostream>
#include <thread>

// // get a rough estimate of how much overhead there is in calling buzy_sleep()
// std::chrono::nanoseconds calc_overhead() {
//     using namespace std::chrono;
//     constexpr size_t tests = 1001;
//     constexpr auto timer = 200;//us

//     auto init = [&timer]() {
//         auto end = steady_clock::now() + timer;
//         while(steady_clock::now() < end);
//     };

//     time_point<steady_clock> start;
//     nanoseconds dur[tests];

//     for(auto& d : dur) {
//         start = steady_clock::now();
//         init();
//         d = steady_clock::now() - start - timer;
//     }
//     std::sort(std::begin(dur), std::end(dur));
//     // get the median value or something a little less as in this example:
//     return dur[tests / 3];
// }

// // initialize the overhead constant that will be used in busy_sleep()
// static const std::chrono::nanoseconds overhead = calc_overhead();

// inline void busy_sleep(std::chrono::nanoseconds t) {
//     auto end = std::chrono::steady_clock::now() + t - overhead;
//     while(std::chrono::steady_clock::now() < end);
// }

int UHD_SAFE_MAIN(int argc, char* argv[])
{

    uhd::device_addr_t my_usrp_ip;
    my_usrp_ip["addr0"] = "192.168.30.2";
    uhd::usrp::multi_usrp::sptr  my_usrp = uhd::usrp::multi_usrp::make(my_usrp_ip);
    std::cout << boost::format("Using Device: %s") % my_usrp->get_pp_string() << std::endl;

    std::cout << "Available GPIO banks: " << std::endl;
    auto banks = my_usrp->get_gpio_banks(0);
    for (auto& bank : banks) {
        std::cout << "* " << bank << std::endl;
    }
    if(false){ // manually set the GPIO pin
        my_usrp->set_gpio_attr("FP0", "CTRL", ATR_CONTROL, ATR_MASKS);
        my_usrp->set_gpio_attr("FP0", "DDR", GPIO_DDR, ATR_MASKS);
        // let's manually set GPIO4 high
        my_usrp->set_gpio_attr("FP0", "OUT", MAN_GPIO_MASK, MAN_GPIO_MASK);
        // finally, let's set up GPIO6 as we described above
        my_usrp->set_gpio_attr("FP0", "ATR_0X", AMP_GPIO_MASK, AMP_GPIO_MASK);
        my_usrp->set_gpio_attr("FP0", "ATR_RX", 0, AMP_GPIO_MASK);
        my_usrp->set_gpio_attr("FP0", "ATR_TX", 0, AMP_GPIO_MASK);//AMP_GPIO_MASK
        // usually, you would want to also make this pin go high when doing
        // full-duplex, but not in this example
        my_usrp->set_gpio_attr("FP0", "ATR_XX", 0, AMP_GPIO_MASK);
    }else{ // switch the 
        // set up some catch-all masks
        uint32_t gpio_line = 0xF; // only the bottom 4 lines: 0xF = 00001111 = Pin 0, 1, 2, 3
        uint32_t all_one = 0xFF;
        uint32_t all_zero = 0x00;

        // reset usrp time to 0.00
        my_usrp->set_time_source("internal");
        my_usrp->set_time_next_pps(uhd::time_spec_t(0.0));
        std::this_thread::sleep_for(std::chrono::milliseconds(2000));

        uhd::time_spec_t now_time = my_usrp->get_time_last_pps();        // define t=0
        // set gpio pins up for output
        my_usrp->set_gpio_attr("FP0", "DDR", all_one, gpio_line, 0);
        my_usrp->set_gpio_attr("FP0", "CTRL", all_zero, gpio_line, 0);
        double Frequency = 42000; //  khz
        double time_step = 1 / Frequency;
        int microseconds_constant = 1000000;
        int sleep_contant_microsec = time_step*microseconds_constant;
        // try{}
        while(true){
            // uhd::time_spec_t now_time = my_usrp->get_time_last_pps();        // define t=0
            for(int i=0;i<Frequency*Frequency;i+=2){
                // set all gpio lines to output 1
                my_usrp->clear_command_time();
                my_usrp->set_command_time(now_time + uhd::time_spec_t(time_step*i));
                my_usrp->set_gpio_attr("FP0", "OUT", all_one, gpio_line, 0);     // reset HIGH (async)
                // set all gpio lines to output 0
                my_usrp->clear_command_time();
                my_usrp->set_command_time(now_time + uhd::time_spec_t(time_step*i+time_step));
                my_usrp->set_gpio_attr("FP0", "OUT", all_zero, gpio_line, 0);    // set LOW @ t=2
                // std::this_thread::sleep_for(std::chrono::microseconds(sleep_contant_microsec));
                // auto end = std::chrono::steady_clock::now() + std::chrono::microseconds(sleep_contant_microsec);
                // while(std::chrono::steady_clock::now() < end);
            }   
            std::cout<<"reach end..."<<std::endl;
            
        }
        
        // // set gpio pins up for output
        // my_usrp->set_gpio_attr("FP0", "OUT", all_one, gpio_line, 0);     // reset HIGH (async)
        // std::cout << "pin HIGH" << std::endl;
        // // set all gpio lines to output 0
        // my_usrp->clear_command_time();
        // my_usrp->set_command_time(now_time + uhd::time_spec_t(2.0));
        // my_usrp->set_gpio_attr("FP0", "OUT", all_zero, gpio_line, 0);    // set LOW @ t=2
        // std::cout << "pin LOW" << std::endl;
        // // set all gpio lines to output 1
        // my_usrp->clear_command_time();
        // my_usrp->set_command_time(now_time + uhd::time_spec_t(4.0));
        // my_usrp->set_gpio_attr("FP0", "OUT", all_one, gpio_line, 0);    // set HIGH @ t=4
        // std::cout << "pin HIGH" << std::endl;
        // // set all gpio lines to output 0
        // my_usrp->clear_command_time();
        // my_usrp->set_command_time(now_time + uhd::time_spec_t(6.0));
        // my_usrp->set_gpio_attr("FP0", "OUT", all_zero, gpio_line, 0);   // set LOW @ t=6
        // std::cout << "pin LOW" << std::endl;
        // my_usrp->clear_command_time();
    }

    return EXIT_SUCCESS;
}

