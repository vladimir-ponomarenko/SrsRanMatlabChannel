#include <iostream>
#include <mutex>
#include <zmq.h>
#include <thread>
#include <cstring>
#include <cmath>
#include <unistd.h>

#include <vector>
#include <random>
#include <complex>

#define BUFFER_MAX 1024 * 1024

static bool run;
std::mutex mutex_send_matlab;

enum class ARGV_CONSOLE {
    ARGV_PORT_1 = 1,
    ARGV_PORT_2,
    ARGV_PORT_1_PROXY,
    ARGV_PORT_2_PROXY,
    ARGV_PORT_API,
    ARGV_MAX
};

void processing_matlab(char *data, int size_data, int size_send, int &size_recv, void *socket_api) {
    std::lock_guard<std::mutex> lock(mutex_send_matlab);

    // printf("Sending %d bytes to MATLAB...\n", size_send);
    int bytes_sent = zmq_send(socket_api, data, size_send, 0);

    if (bytes_sent == -1) {
        perror("!!! zmq_send to MATLAB failed");
        fprintf(stderr, "!!! Error sending to MATLAB: %s\n", zmq_strerror(zmq_errno()));
        size_recv = -1;
        return;
    }
    if (bytes_sent != size_send) {
         fprintf(stderr, "!!! Warning: zmq_send sent %d bytes instead of %d\n", bytes_sent, size_send);
    }

    // printf("Waiting for reply from MATLAB...\n");
    size_recv = zmq_recv(socket_api, data, size_data, 0);

    if (size_recv == -1) {
        perror("!!! zmq_recv from MATLAB failed");
        int err = zmq_errno();
        fprintf(stderr, "!!! Error receiving from MATLAB: %s (%d)\n", zmq_strerror(err), err);
        if (err == EAGAIN) {
            fprintf(stderr, "!!! MATLAB did not respond within timeout.\n");
        }
    } else {
         // printf("Received %d bytes from MATLAB.\n", size_recv); 
    }
}

int iter = 0;


void processing_data(char *data, int size_data, int size_send, int &size_recv, void *socket_api) {

    processing_matlab(data, size_data, size_send, size_recv, socket_api);
    // ...
}

void thread_proxy(void *zrecv, void *zsend, void *socket_api, int id) {
    
    char buffer[BUFFER_MAX];
    int size;
    int size_recv;
    printf("[%d] start\n", id);
    while(1) {
        memset(buffer, 0, sizeof(buffer));
        size = zmq_recv(zrecv, buffer, sizeof(buffer), 0);
        
        if(size == -1) {
            continue;
        }
        if(size > 1000)
            processing_data(buffer, BUFFER_MAX, size, size_recv, socket_api);
        
        zmq_send(zsend, buffer, size, 0);
    }
}

void thread_proxy_2(void *zrecv, void *zsend, void *socket_api, int id) {
    
    char buffer[BUFFER_MAX];
    int size;
    int size_recv;
    printf("[%d] start\n", id);
    while(1) {
        memset(buffer, 0, sizeof(buffer));
        size = zmq_recv(zsend, buffer, sizeof(buffer), 0);
        if(size == -1) {
            continue;
        }
        if(size > 1000)
            processing_data(buffer, BUFFER_MAX, size, size_recv, socket_api);
        zmq_send(zrecv, buffer, size, 0);
    }
}


int main(int argc, char *argv[]){

    if(argc < static_cast<int>(ARGV_CONSOLE::ARGV_MAX) - 1) {
        printf("Error: not found argv\n");
        return -1;
    }
    {
        FILE *file = fopen("../statistics/statistics.txt", "w");
        if(file) {
            fclose(file);
        }
    }
    int port1 = std::stoi(argv[static_cast<int>(ARGV_CONSOLE::ARGV_PORT_1)]);
    int port2 = std::stoi(argv[static_cast<int>(ARGV_CONSOLE::ARGV_PORT_2)]);
    int port1_proxy = std::stoi(argv[static_cast<int>(ARGV_CONSOLE::ARGV_PORT_1_PROXY)]);
    int port2_proxy = std::stoi(argv[static_cast<int>(ARGV_CONSOLE::ARGV_PORT_2_PROXY)]);
    int port_api = std::stoi(argv[static_cast<int>(ARGV_CONSOLE::ARGV_PORT_API)]);
    
    void *context = zmq_ctx_new ();
    void *requester = zmq_socket (context, ZMQ_REQ);
    void *requester2 = zmq_socket (context, ZMQ_REQ);
    void *requester_api = zmq_socket (context, ZMQ_REQ);
    printf("%d %d %d %d\n", port1, port2, port1_proxy, port2_proxy);
    std::string addr_recv_1 = "tcp://localhost:" + std::to_string(port1);
    std::string addr_recv_2 = "tcp://localhost:" + std::to_string(port2);
    std::string addr_send_1 = "tcp://*:" + std::to_string(port1_proxy);
    std::string addr_send_2 = "tcp://*:" + std::to_string(port2_proxy);
    std::string addr_socket_api = "tcp://localhost:" + std::to_string(port_api);
    printf("connecting to %s...\n", addr_recv_1.c_str());
    zmq_connect (requester, addr_recv_1.c_str());
    printf("connecting to %s...\n", addr_recv_2.c_str());
    zmq_connect (requester2, addr_recv_2.c_str());
    printf("connecting to %s...\n", addr_socket_api.c_str());
    zmq_connect (requester_api, addr_socket_api.c_str());
    void *socket_proxy1 = zmq_socket(context, ZMQ_REP);
    void *socket_proxy2 = zmq_socket(context, ZMQ_REP);
    if(!socket_proxy1) {
        perror("zmq_socket");
        return 1;
    }
    if(!socket_proxy2) {
        perror("zmq_socket");
        return 1;
    }
    
    if(zmq_bind(socket_proxy1, addr_send_1.c_str()) != 0) {
        perror("zmq_bind");
        return 1;
    }
    if(zmq_bind(socket_proxy2, addr_send_2.c_str())) {
        perror("zmq_bind");
        return 1;
    }
    printf("[Init]\n");
    std::thread thr1;
    std::thread thr2;
    std::thread thr3;
    std::thread thr4;
    thr1 = std::thread(thread_proxy, socket_proxy1, requester, requester_api, 1);
    thr1.detach();
    thr3 = std::thread(thread_proxy_2, socket_proxy1, requester, requester_api, 2);
    thr3.detach();
    thr2 = std::thread(thread_proxy, socket_proxy2, requester2, requester_api, 3);
    thr2.detach();
    thr4 = std::thread(thread_proxy_2, socket_proxy2, requester2, requester_api, 4);
    thr4.detach();
    while(1) {
        sleep(1);
    }
    printf("End client\n");
    zmq_close (requester);
    zmq_close (requester2);
    zmq_close (socket_proxy1);
    zmq_close (socket_proxy2);
    zmq_close (requester_api);
    zmq_ctx_destroy (context);
    printf("[Clear]\n");
    return 0;
}