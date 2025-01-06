import websockets, asyncio, json


async def connect_and_send(rule_name, status="update"):
    uri = "ws://localhost:8009/xplorecov/wsa/backend/"

    try:
        async with websockets.connect(uri) as websocket:
            print(f"Connected to WebSocket server - {status}")

            response = await websocket.recv()
            print(f"Received: {response}")

            message = {"type": rule_name, "status": status}
            await websocket.send(json.dumps(message))
            print(f"Sent status: {status}")

            # Close connection
            await websocket.close()
            print("Connection closed")

    except Exception as e:
        print(f"WebSocket error: {str(e)}")


def run_websocket_message(rule_name, status="update"):
    """Synchronous wrapper for async websocket connection"""
    # Create new event loop for this thread
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    try:
        loop.run_until_complete(connect_and_send(rule_name, status))
    finally:
        loop.close()


rule all:
    input:
        # f"{config['OutputDir']}/result/nextclade/clade_report.tsv",
        # f"{config['OutputDir']}/result/pangolin-usher/lineage.tsv",
        f'{config["OutputDir"]}/result/combined_report.tsv',


include: "rules/annotation/nextclade.smk"
include: "rules/annotation/pangolin-usher.smk"
include: "rules/combine/index.smk"


onstart:
    run_websocket_message("WORKFLOW", "Started Workflow")


onsuccess:
    run_websocket_message("WORKFLOW", "Completed Workflow")


onerror:
    run_websocket_message("WORKFLOW", "Cancelled Workflow")
