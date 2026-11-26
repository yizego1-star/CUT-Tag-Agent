"""A bounded observe/decide/act/re-observe loop over two real tool rounds."""
import copy
from datetime import datetime, timezone

from .model import ModelError, minimal_summary, validate_decision


def drive_agent(run_id, sample_id, coverage, model, cancelled, update):
    rounds, decisions, requests = [], [], []
    size = 25

    def snapshot():
        return {"rounds": copy.deepcopy(rounds), "decisions": copy.deepcopy(decisions), "model_requests": copy.deepcopy(requests),
                "actual_tool_executions": sum(r.get("actual_tool_executions", 0) for r in rounds),
                "actual_model_requests": sum(bool(r.get("http_attempted")) for r in requests),
                "reruns_completed": int(len(rounds) == 2 and rounds[-1]["status"] == "completed")}

    def end(status, reason, **extra):
        return {**snapshot(), "status": status, "end_reason": reason, **extra}

    for index in range(2):
        if cancelled():
            return end("cancelled", "用户停止任务，未继续调用工具或模型")
        key = "baseline" if index == 0 else "rerun"
        update(key, "正在运行真实基线计算" if index == 0 else "正在执行模型建议的一次真实重跑", snapshot())
        try:
            output = coverage.run(run_id, sample_id, cancelled, lambda *_: None, bin_size=size, round_index=index)
        except Exception:
            return end("failed", "真实工具运行未完成；请查看本地工具环境与日志")
        rounds.append({"index": index, "parameters": {"bin_size": size}, **output})
        update(key, "工具轮次已返回", snapshot())
        if output["status"] != "completed":
            return end(output["status"], output.get("execution", {}).get("error") or "工具未完成，未继续请求模型")
        if cancelled():
            return end("cancelled", "用户停止任务，未发送下一次模型请求")
        summary = minimal_summary(output["result"], output["bam_stats"], index)
        request = {"round_index": index, "summary": summary, "http_attempted": False,
                   "created_at": datetime.now(timezone.utc).isoformat()}
        requests.append(request)
        key = "decision" if index == 0 else "review"
        update(key, "模型正在读取获准的最小摘要并提出建议" if index == 0 else "模型正在读取重跑摘要并作结束判断", snapshot())

        def sent():
            if cancelled():
                raise ModelError("任务已停止，未发送模型请求")
            request["http_attempted"] = True
            update(key, "已向模型发送获准的最小摘要", snapshot())

        try:
            answer = model.decide(summary, on_request=sent)
            decision = validate_decision(answer["decision"], size, 1 - index)
        except ModelError as exc:
            request["error"] = str(exc)
            return end("cancelled" if cancelled() else "failed", str(exc), error=str(exc))
        request["usage"] = answer.get("usage", {})
        decisions.append({"round_index": index, "model": answer["model"], **decision})
        update(key, "模型建议已解析并通过动作范围检查", snapshot())
        if cancelled():
            return end("cancelled", "任务已停止，在途模型回复未触发后续计算")
        if decision["action"] == "finish":
            return end("completed_agent", decision["reason"])
        size = decision["bin_size"]
    return end("failed", "已达到一次重跑上限，未继续执行")
