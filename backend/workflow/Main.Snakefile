from datetime import datetime
import websockets, asyncio, json, os, pytz


async def connect_and_send(step_name, step_id, status, message_type):
    timezone = pytz.timezone("Asia/Kolkata")
    uri = f"ws://10.10.6.80/xplorecov/ws/analysis/{config['UserID']}/{config['AnalysisID']}/?BACKEND_WEBSOCKET_UUID={os.getenv('BACKEND_WEBSOCKET_UUID')}"
    try:
        async with websockets.connect(uri) as websocket:
            response = await websocket.recv()
            message = {
                "message_type": message_type,
                "message": (
                    {
                        "status": status,
                        "step_id": step_id,
                        "step_name": step_name,
                        "timestamp": datetime.now(timezone).isoformat(),
                    }
                    if message_type == "analysis_update"
                    else {
                        "status": status,
                        "timestamp": datetime.now(timezone).isoformat(),
                    }
                ),
            }
            await websocket.send(json.dumps(message))
            await websocket.close()
    except Exception as e:
        print(f"WebSocket error: {str(e)}")


def run_websocket_message(step_name, step_id, status, message_type="analysis_update"):
    """Synchronous wrapper for async websocket connection"""
    # Create new event loop for this thread
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    try:
        loop.run_until_complete(
            connect_and_send(step_name, step_id, status, message_type)
        )
    finally:
        loop.close()


rule all:
    input:
        f'{config["OutputDir"]}/result/combined_report.tsv',
        f'{config["OutputDir"]}/result/combined_report.feather',


include: "rules/annotation/nextclade.smk"
include: "rules/annotation/pangolin-usher.smk"
include: "rules/combine/index.smk"


onstart:
    run_websocket_message("Queuing Analysis", "step3", "end")


onsuccess:
    # SEND COMPLETE TO OVERALL STATUS
    run_websocket_message(
        "WORKFLOW", "Completed Workflow", "SUCCESS", "workflow_update"
    )


onerror:
    # SEND ERROR TO OVERALL STATUS
    run_websocket_message("WORKFLOW", "Error Workflow", "ERROR", "workflow_update")
